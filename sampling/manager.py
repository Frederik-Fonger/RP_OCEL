import csv
import math
import os
import time
from copy import deepcopy, copy
from datetime import datetime
import ntpath

import pandas as pd
import pm4py
from fontTools.ttLib.woff2 import base128Size

import sampling.base_algorithms_AB_RP

from ocpa.objects.log.importer.ocel2.xml import factory as ocel_import_factory_XML
from ocpa.objects.log.importer.ocel import factory as ocel_import_factory_JSON

from sampling import base_algorithms_AB_RP
from scipy.optimize import bisect, minimize_scalar

"""
OCEL_sampling

Class for one sampling run. With the .apply() function the sampling started. 
"""





class OCEL_sampling:


    def __init__(self):
        """
        Stores the variable need while sampling. The variable are set, as soon as there are available.
        """
        self.ocel_pm4py = None
        self.ocel_pm4py_orginal = None
        self.ocel_ocpa = None
        self.ocel_ocpa_orginal = None
        self.connectivity_degree = None
        self.number_of_object_conections = None
        self.sampling_algo = None
        self.sampling_ratio = None
        self.connectivity_threshold = None
        self.pe_list_DFRs = None
        self.event_count_origial_ocel = None
        self.target_event_count_in_sample = None
        self.list_of_selected_objects = None
        self.path_for_result = None



    def import_ocel_xml(self, path, file_type):
        """
        Function for importing ocel XML files and Storing them is the OCEL_sampling Object.
        :param path: file path of OCEL-XML
        :param file_type: "XML"

        :return: None
        """
        if file_type == "XML":
            self.ocel_pm4py = pm4py.read.read_ocel2_xml(path)
            self.ocel_pm4py_orginal = pm4py.read.read_ocel2_xml(path)
            self.ocel_ocpa = ocel_import_factory_XML.apply(path)
            self.ocel_ocpa_orginal = ocel_import_factory_XML.apply(path)
            self.event_count_origial_ocel = len(self.ocel_pm4py.events)

        else:
            raise TypeError("File type " + file_type + " not supported")
        print("import of " + path + " successful")

    def calculate_connectivity_per_ot(self):
        """
        This function calculates the average number of related object for each object type. Two object are only related
        once, even if they have interaction in multiple events.
        :return:
        connectivity_degree:            A dataframe with the number of average object relations for each object type
        number_of_object_conections:    The total number of object relation for an object type
        """
        # import the OCEL


        # get the information table about the object interactions provides by pm4py
        object_summary = pm4py.ocel.ocel_objects_interactions_summary(self.ocel_pm4py)
        # combine the two object ids from an interaction to generate a uniqe identifier for the relation of this two objects
        object_summary['combined_oids'] = object_summary['ocel:oid'] + "_" + object_summary['ocel:oid_2']
        # delete all dublicates of the uniqe identifier for the relation of this two objects, as each object to object
        # realated should be counted once. The listed for each event with the two activites again.
        object_summary = object_summary.drop_duplicates(subset=['combined_oids'], keep='first')
        # delete rows the object-to-object relation is between objects of the same type
        # object_summary = object_summary[object_summary['ocel:type'] != object_summary['ocel:type_2']]
        # count the number ob object relations per object type
        number_of_object_conections = object_summary['ocel:type'].value_counts()
        print("number objects connections:")
        print(number_of_object_conections)
        # count the number ob objects per type
        number_objects_per_type = self.ocel_pm4py.objects['ocel:type'].value_counts()
        print("number objects per type:")
        print(number_objects_per_type)

        connectivity_degree = []
        # calculate for each object type the average number of object to object relations
        for object_type in number_of_object_conections.index:
            connectivity_of_one_object = number_of_object_conections[object_type] / number_objects_per_type[object_type]
            connectivity_degree.append([object_type, connectivity_of_one_object])

        # make a pandas df from the before used list and name the columns correctly
        connectivity_degree = pd.DataFrame(connectivity_degree)
        connectivity_degree = connectivity_degree.rename(columns={0: 'object_type', 1: 'average_related_objects'})
        print(connectivity_degree)

        self.connectivity_degree = connectivity_degree
        self.number_of_object_conections = number_of_object_conections
        print("Connectivity per object type calculated")
        return None

    def calculate_connectivity_score_v1(self):
        """
        This function calculates a connectivity score for each object based on the number of connections from object to
        object. This version uses a log 2 and divides the connection from each object by the number of connection from
        the max object. :param connectivity_degree: :param number_of_object_conections: :return:
        """

        def log2_or_zero(x):
            print(x)
            if x <= 1:
                return 0
            else:
                return math.log2(x)
        # apply the log2 in the number of average object relation for each object type
        self.connectivity_degree['average_related_objects_log'] = self.connectivity_degree['average_related_objects'].apply(
            log2_or_zero)

        # get the highest number of average object relation
        max_count_of_average_related_objects = self.connectivity_degree['average_related_objects_log'].max()
        # divide the number of object relation for each objeyt type by the highest number of average object relation over all object types
        self.connectivity_degree['score'] = self.connectivity_degree[
                                           'average_related_objects_log'] / max_count_of_average_related_objects
        print("Connectivity score V1 is calculated")


    def filter_object_types_by_connectivity_threshold(self):
        """
        the pm4py-OCEL in the OCELSampling object is filtered by the object types that are below the set connectivity
        threshold. The ocel_pm4py is updated and the OCPA-ocel is also updated. :return: None
        """
        print("Connectivity degree of the object types: ")
        print(self.connectivity_degree[['object_type', 'score']])

        # create a list of all objects of that belog to the object type that are below the set connectivity
        # threshold.
        self.list_of_selected_objects = self.connectivity_degree['object_type'].loc[self.connectivity_degree['score'] <= self.connectivity_threshold].to_list()
        self.ocel_pm4py = pm4py.filtering.filter_ocel_object_types(self.ocel_pm4py, self.list_of_selected_objects)
        OCEL_sampling.update_from_pm4py_to_ocpa(self)
        print("OCEL has been filtered by connectivity threshold")


    def calculate_DFRs_per_process_execution(self):
        """
        calculate the object-centic DFRs. THis is not the direct successor in the event log, but the successor for
        the event of the same object.
        :return: None
        """

        pe_list_DRFs = []

        # iterate though the PEs (pe is a list of events)
        for pe in self.ocel_ocpa.process_executions:

            # the PE is created in a pm4py ocel, by filtering all event out that not belong to the PE
            filtered_ocel = pm4py.filtering.filter_ocel_events(self.ocel_pm4py, list(pe))
            DFRs = {}

            # to the DFRs for each object individually, the longs need to be flattened.
            for object_type in pm4py.ocel.ocel_get_object_types(filtered_ocel):

                # calculatin the DFRs in the flat log
                flat_log = pm4py.ocel.ocel_flattening(filtered_ocel, object_type)
                dfg = pm4py.discovery.discover_dfg(flat_log)
                DFRs = {k: DFRs.get(k, 0) + dfg[0].get(k, 0) for k in set(DFRs) | set(dfg[0])}

            pe_list_DRFs.append(DFRs)
        self.pe_list_DFRs = pe_list_DRFs
        print("DFRs per process execution are calculated")

    def sampling_and_create_list_of_selected_eventIDs(self):
        """
        function to select the chosen sampling algorithm.

        :return: a list of chosen process execution for the sample
        """

        if self.sampling_algo == "AB":
            return sampling.base_algorithms_AB_RP.ocel_allbehaviour_sampling(self)
        elif self.sampling_algo == "RP":
            return sampling.base_algorithms_AB_RP.ocel_remainder_plus_sampling(self)
        elif self.sampling_algo == "None":
            return list(range(len(self.pe_list_DFRs)))
        else:
            raise Exception("sampling method: " + self.sampling_algo + " not supported")


    def update_from_pm4py_to_ocpa(self):
        """
        This function updates changes from the pm4py-ocel object to the ocpa-ocel object. This is done by save and load of a temp XML file.
        :return:
        """
        pm4py.write.write_ocel2_xml(self.ocel_pm4py, "temp.xmlocel")
        self.ocel_ocpa = ocel_import_factory_XML.apply("temp.xmlocel")

    def update_from_ocpa_to_pm4py(self):
        # Todo this function
        pass

    def export_pm4py_ocel(self, pm4py_ocel, file_name, create_folder_with_sample_and_meta_data, file_path, original_sampling_ratio):
        if not create_folder_with_sample_and_meta_data:
            pm4py.write.write_ocel2_xml(pm4py_ocel, "sample" + ".xmlocel")
        else:
            # Get the current time as a string
            folder_name = 'output/' +datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + "_ocel_sampling"

            # Create the folder
            os.makedirs(folder_name)

            file_name = ntpath.basename(file_path)
            file_name = os.path.splitext(file_name)[0]
            self.path_for_result = folder_name + '/'
            pm4py.write.write_ocel2_xml(pm4py_ocel, self.path_for_result + file_name + "_sample_" + self.sampling_algo + "_" + str(original_sampling_ratio) + ".xmlocel")


    def save_and_print_meta_data(self, create_folder_with_sample_and_meta_data, sampling_time, path, original_sampling_ratio):
        print("Object types used for selecting events: ")
        print(self.list_of_selected_objects)
        if create_folder_with_sample_and_meta_data:
            self.connectivity_degree['used_for_selecting'] = self.connectivity_degree['object_type'].isin(
                self.list_of_selected_objects)

            self.connectivity_degree.to_csv(self.path_for_result + 'object_types.csv', sep=';')

            meta_data = {}
            meta_data['connectivity_threshold'] = self.connectivity_threshold
            meta_data['sampling_time'] = sampling_time
            meta_data['object types for selection'] = self.list_of_selected_objects
            meta_data['file path'] = path
            meta_data['sampling algo'] = original_sampling_ratio
            meta_data['sample ratio for step one'] = self.sampling_algo
            meta_data['sampling ratio'] = self.sampling_ratio


            with open(self.path_for_result + 'meta_data.csv', mode="w", newline="") as file:
                w = csv.DictWriter(file, meta_data.keys(), delimiter=";")
                w.writeheader()
                w.writerow(meta_data)


    def create_sample_from_list_of_PEs(self, remainder_plus_sample):
        """
        This function creating a sample in pm4py ocel format from a list of process execution.

        :param remainder_plus_sample: list of process execution from the sample

        :return: sample in pm4py ocel format
        """

        # annotate the process execution to every event in the log, so that it is possible to filter the log by process executions.
        self.ocel_ocpa.log.log['process_execution'] = self.ocel_ocpa.log.log.index.map(
            self.ocel_ocpa.process_execution_mappings)

        # the number of the process execution is not an int but a list with one int inside. This funktion anpacks the list
        # print(self.ocel_ocpa_orginal.log.log['process_execution'])
        self.ocel_ocpa.log.log['process_execution'] = self.ocel_ocpa.log.log[
            'process_execution'].apply(lambda x: x[0])

        # create a ocel sample (dataframe from inside the ocpa-obbject) with only the chosen PEs
        sample_ocpa_dataframe = self.ocel_ocpa.log.log[self.ocel_ocpa.log.log['process_execution'].isin(remainder_plus_sample)]

        # to convert from ocpa format to pm4py format, a list of selected eventIDs in created.
        list_of_selected_event_ids_for_sample = list(sample_ocpa_dataframe['event_id'])

        # cerate a sample in pm4py ocel format by filtering the pm4py ocel with the selected eventIDs
        filtered_ocel = pm4py.filtering.filter_ocel_events(self.ocel_pm4py_orginal,
                                                           list_of_selected_event_ids_for_sample)


        # FROM HERE NEW FOR ADDING RELATED EVENT OUTSIDE OF PEs (if not all event types are in sample) -----------------------------------------

        def check_if_all_event_types_are_in_sample(self, filtered_ocel):
            sample_event_types = filtered_ocel.events['ocel:activity'].unique()
            original_event_types = self.ocel_pm4py_orginal.events['ocel:activity'].unique()
            test = set(sample_event_types) == set(original_event_types)
            return set(sample_event_types) == set(original_event_types)

        if not check_if_all_event_types_are_in_sample(self, filtered_ocel):


            objects_in_sample_from_not_used_object_types = filtered_ocel.objects[~filtered_ocel.objects['ocel:type'].isin(self.list_of_selected_objects)]

            # get list from object that are not in the sample
            list_of_not_used_objects = self.ocel_pm4py_orginal.objects.merge(filtered_ocel.objects, how='outer', indicator=True).query('_merge == "left_only"').drop(columns=['_merge'])
            # get the object that are in the sample and not in the cut away part of the ocel, but are not from the object types that are selected for sampling
            list_of_objects_only_in_sample_of_not_used_type = objects_in_sample_from_not_used_object_types.merge(list_of_not_used_objects, how='outer', indicator=True).query('_merge == "left_only"').drop(columns=['_merge'])

            # get the event ids from the events that are related to the to the object that are in the sample and not in the cut away part of the ocel, but are not from the object types that are selected for sampling
            event_ids_from_list_of_objects_only_in_sample_of_not_used_type = self.ocel_pm4py_orginal.relations[self.ocel_pm4py_orginal.relations['ocel:oid'].isin(list_of_objects_only_in_sample_of_not_used_type['ocel:oid'])]
            event_ids_from_list_of_objects_only_in_sample_of_not_used_type = event_ids_from_list_of_objects_only_in_sample_of_not_used_type.drop_duplicates(subset=['ocel:eid'])
            # temp_list = event_ids_from_list_of_objects_only_in_sample_of_not_used_type['ocel:eid'].unique().tolist()
            list_of_selected_event_ids_for_sample = list_of_selected_event_ids_for_sample + event_ids_from_list_of_objects_only_in_sample_of_not_used_type['ocel:eid'].unique().tolist()

            list_of_selected_event_ids_for_sample = list(set(list_of_selected_event_ids_for_sample))

            # to convert from ocpa format to pm4py format, a list of selected eventIDs in created.
            # list_of_selected_event_ids_for_sample = list(sample_ocpa_dataframe['event_id'])

            # cerate a sample in pm4py ocel format by filtering the pm4py ocel with the selected eventIDs
            filtered_ocel = pm4py.filtering.filter_ocel_events(self.ocel_pm4py_orginal,
                                                               list_of_selected_event_ids_for_sample)

        # END OF NEW ADDITION -------------------------------------

        return filtered_ocel


    def striped_sampling(self, ocel_pm4py, ocel_ocpa, sample_ratio, sampling_algo="AB", connectivity_threshold=0.5):

        self.sampling_ratio = sample_ratio
        self.sampling_algo = sampling_algo
        self.connectivity_threshold = connectivity_threshold

        self.ocel_pm4py = ocel_pm4py
        self.ocel_pm4py_orginal = deepcopy(ocel_pm4py)
        self.ocel_ocpa = ocel_ocpa
        self.ocel_ocpa_orginal = deepcopy(ocel_ocpa)
        self.event_count_origial_ocel = len(self.ocel_pm4py.events)

        # Calculate the start time
        start = time.time()

        # Connectivity per object type
        OCEL_sampling.calculate_connectivity_per_ot(self)
        OCEL_sampling.calculate_connectivity_score_v1(self)

        # create a OCEl with only object-type below threshold
        OCEL_sampling.filter_object_types_by_connectivity_threshold(self)

        # calcualte DFRs per PE
        OCEL_sampling.calculate_DFRs_per_process_execution(self)

        # call RP or AB with PE and DFRs (return eventIDs)
        list_of_selected_PEs = OCEL_sampling.sampling_and_create_list_of_selected_eventIDs(self)

        # construct the sample from the list of chosen process executions.
        sample_ocel_pm4py = OCEL_sampling.create_sample_from_list_of_PEs(self, list_of_selected_PEs)

        end = time.time()
        sampling_time = end - start

        print("sampling finished successfully!")

        return sample_ocel_pm4py


    def apply(self, file_path, sample_ratio, file_type= "XML", sampling_algo="AB", connectivity_threshold=0.5, create_folder_with_sample_and_meta_data = False):
        """
        Function the start the OCEL sampling.

        :param file_path: path to the OCEL file
        :param sample_ratio: selected sample ratio for sample (e.g. 0.5 for 50% of the events from original log)
        :param file_type: e.g. "XML
        :param sampling_algo: "RP" for Remainder Plus or "AB" for AllBehaviour
        :param connectivity_threshold: A threshold with objects are used for selection the events. Object with a higher
                    connectivity degree then the threshold are no take into account.
        :return: None
        """

        self.sampling_ratio = sample_ratio
        self.sampling_algo = sampling_algo
        self.connectivity_threshold = connectivity_threshold

        # import file
        OCEL_sampling.import_ocel_xml(self, file_path, file_type)
        self.target_event_count_in_sample = copy(int(self.event_count_origial_ocel * self.sampling_ratio))
        # Calculate the start time
        start = time.time()

        # Connectivity per object type
        OCEL_sampling.calculate_connectivity_per_ot(self)
        OCEL_sampling.calculate_connectivity_score_v1(self)

        # create a OCEl with only object-type below threshold
        OCEL_sampling.filter_object_types_by_connectivity_threshold(self)

        # calcualte DFRs per PE
        OCEL_sampling.calculate_DFRs_per_process_execution(self)
        
        #-----------------Start of sample size search ----------------------
        target_event_count_in_sample = copy(self.target_event_count_in_sample)
        
        # inizialize values that are used in search function
        sample_ocel_pm4py = None
        list_of_selected_PEs = []
        
        def search_fitting_sampleratio(sample_ratio):
            """
            Function to search for the best fitting sample ratio. This must be done because the first step of the sampling algorithm can take a sample ratio, but the second step adds an unknonw number of events. To find the fitting sample ratio for the first step, the event count in the end is used to optimize towards a fitting initial sample ratio.
            :param sample_ratio: sample ratio for this try
            :return: the absolute difference between the target event count and the current event count in the sample.
            """
            # set the sample ratio for the next sampling run 
            self.sampling_ratio = sample_ratio
            # call RP or AB with PE and DFRs (return eventIDs)
            list_of_selected_PEs = OCEL_sampling.sampling_and_create_list_of_selected_eventIDs(self)
            # construct the sample from the list of chosen process executions.
            sample_ocel_pm4py = OCEL_sampling.create_sample_from_list_of_PEs(self, list_of_selected_PEs)
            # returen the absolute difference between the target event count and the current event count in the sample.
            return  abs(len(sample_ocel_pm4py.events) - target_event_count_in_sample)

        # track time for optimization
        start_opi = time.time()
        # find a sample ratio for the initial step, so that after the second step the events count corresponds to the origial sample ratio.
        ideal_sampleratio = minimize_scalar(search_fitting_sampleratio, bounds=(0,sample_ratio + 0.05),method='bounded' ).x  # bounds you know are safe
        end_opi = time.time()

        self.sampling_ratio = ideal_sampleratio
        # list_of_selected_PEs = OCEL_sampling.sampling_and_create_list_of_selected_eventIDs(self)
        # sample_ocel_pm4py = OCEL_sampling.create_sample_from_list_of_PEs(self, list_of_selected_PEs)

        # -----------------end of sample size search ----------------

        # call RP or AB with PE and DFRs (return eventIDs)
        list_of_selected_PEs = OCEL_sampling.sampling_and_create_list_of_selected_eventIDs(self)
        # construct the sample from the list of chosen process executions.
        sample_ocel_pm4py = OCEL_sampling.create_sample_from_list_of_PEs(self, list_of_selected_PEs)

        end = time.time()
        sampling_time = end - start

        OCEL_sampling.export_pm4py_ocel(self, sample_ocel_pm4py, "result", create_folder_with_sample_and_meta_data, file_path, sample_ratio)

        OCEL_sampling.save_and_print_meta_data(self,create_folder_with_sample_and_meta_data, sampling_time, file_path, sample_ratio)
        print("Time for optimization: ", end_opi - start_opi)
        print("Sampling runtime: " + str(sampling_time))
        print("Event count: " + str(len(sample_ocel_pm4py.events.index)))
        print("finished successfully!")


