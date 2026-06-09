import time

import pandas as pd
import pm4py
from collections import Counter
import os

# from lxml.parser import result

import sampling.manager as smpl
from ocpa.objects.log.importer.ocel2.xml import factory as ocel_import_factory
from ocpa.algo.discovery.ocpn import algorithm as ocpn_discovery_factory
from ocpa.algo.conformance.precision_and_fitness import evaluator as quality_measure_factory
# from ocpa.algo.enhancement.token_replay_based_performance import algorithm as token_based_replay_algorithm
from ocpa_main.ocpa.algo.conformance.token_based_replay import algorithm as token_based_replay_algorithm
from ocpa.visualization.oc_petri_net import factory as ocpn_vis_factory

class Evaluation:


    def __init__(self):
        """
        Stores the variable need while evaluating. The variable are set, as soon as there are available.
        """
        self.sample_pm4py = None
        self.sample_ocpa = None
        self.ocel_pm4py = None
        self.ocel_pm4py_orginal = None
        self.ocel_ocpa = None
        self.ocel_ocpa_orginal = None
        self.sampling_algo = None
        self.sampling_ratio = None
        self.quality_metrics = None
        self.token_replay = False



    def import_ocel(self, path):
        self.ocel_pm4py = pm4py.read.read_ocel2_xml(path)
        self.ocel_pm4py_orginal = pm4py.read.read_ocel2_xml(path)
        self.ocel_ocpa = ocel_import_factory.apply(path)
        self.ocel_ocpa_orginal = ocel_import_factory.apply(path)


    def update_sample_from_pm4py_to_ocpa(self):
        """
        This function updates changes from the pm4py-ocel object to the ocpa-ocel object. This is done by save and load of a temp XML file.
        :return:
        """
        pm4py.write.write_ocel2_xml(self.sample_pm4py, "temp.xmlocel")
        self.sample_ocpa = ocel_import_factory.apply("temp.xmlocel")

    def update_from_ocpa_to_pm4py(self):
        # Todo this function
        pass

    def calculate_quality_metrics(self, ocel, ocpn):

        if self.token_replay:
            metrics_dict = token_based_replay_algorithm.apply(ocel, ocpn)
        else:
            precision, fitness = quality_measure_factory.apply(ocel, ocpn)
            metrics_dict = {
                'fitness' : fitness,
                'precision': precision
            }

        self.quality_metrics = metrics_dict


    def calculate_DFRs(self, ocap_log):
        """
        calculate the object-centic DFRs. THis is not the direct successor in the event log, but the successor for
        the event of the same object.
        :return: None
        """

        pe_list_DRFs = []

        # iterate though the PEs (pe is a list of events)
        for pe in ocap_log.process_executions:

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

        # Sum the PEs up
        DFRs_of_OCEL = Counter()
        for pe in pe_list_DRFs:
            DFRs_of_OCEL.update(pe)

        # Convert back to a regular dict if needed
        DFRs_of_OCEL = dict(DFRs_of_OCEL)

        return DFRs_of_OCEL

    def calculate_MAE_and_coverage(self, sample, original_log):

        # import log
        Evaluation.import_ocel(self, original_log)
        self.sample_ocpa = ocel_import_factory.apply(sample)

        dfrs_sample = self.calculate_DFRs(self.sample_ocpa)
        dfrs_original = self.calculate_DFRs(self.ocel_ocpa_orginal)
        print("test")
        difference = {k: v for k, v in dfrs_sample.items() if k not in dfrs_original}
        1
        difference.update({k: v for k, v in dfrs_original.items() if k not in dfrs_sample})

        dfrs_num_of_dfrs_sample_that_are_also_in_original = len(dfrs_sample) - len(difference)

        sum_sample = sum(dfrs_sample.values())
        sum_original = sum(dfrs_original.values())
        sum_difference = sum(difference.values())

        print("DFRs sampled: " + str(dfrs_sample))
        print("DFRs original: " + str(dfrs_original))
        print("difference: " + str(difference))

        print("Number of DFRs in sample: " + str(len(dfrs_sample)))
        print("Number of DFRs in original: " + str(len(dfrs_original)))
        print("Number of DFRs that are also in original: " + str(dfrs_num_of_dfrs_sample_that_are_also_in_original))

        coverage = dfrs_num_of_dfrs_sample_that_are_also_in_original / len(dfrs_original)

        # MAE
        # calcualte the unification of all keys
        all_keys = dfrs_sample.keys() | dfrs_original.keys()  # union of keys

        # calculate the error (difference) of each DFR und use 0 as a frequency if the DFR dose not exist in the log
        errors = [abs(dfrs_sample.get(k, 0) - dfrs_original.get(k, 0)) for k in all_keys]
        # divide the sum of the error by the number of idividula DFRs to calculate the MAE
        mae = sum(errors) / len(all_keys)

        #
        changed_some_column, total_in_both, totally_unchanged, percentage_of_changed_columns = Evaluation.calculate_obj_difference(self)

        # mae =
        results = {"Sample": sample, "original log": original_log, "MAE": mae, "coverage": coverage, "changed_some_column": changed_some_column, "total_in_both": total_in_both, "totally_unchanged": totally_unchanged, "percentage_of_changed_columns": percentage_of_changed_columns}
        print("MAE: " + str(mae))
        print("coverage: " + str(coverage))
        print("changed_some_column: " + str(changed_some_column))
        print("total_in_both: " + str(total_in_both))
        print("totally_unchanged: " + str(totally_unchanged))
        print("percentage_of_changed_columns: " + str(percentage_of_changed_columns))
        print(results)

        pd.DataFrame.from_dict(results, orient='index', columns=['value']).to_csv("MAE_and_coverage_"+ sample.rsplit("/", 1)[1] +".csv")

        return coverage, mae

    def calculate_obj_difference(self):

        # assume df_full is the larger DataFrame, df_sub is the subset-with-changes
        # both have an "eventID" column

        # 1. index both by eventID so we can align rows easily
        df1 = self.ocel_ocpa.log.log.set_index('event_id')
        df2 = self.sample_ocpa.log.log.set_index('event_id')

        # 2. find the eventIDs present in both
        common_ids = df1.index.intersection(df2.index)

        # 3. extract only those common rows
        df1_common = df1.loc[common_ids]
        df2_common = df2.loc[common_ids]

        # 4. compare the rows element‑wise
        #    this gives a DataFrame of booleans: True where values match
        matches = df1_common.eq(df2_common)

        # 5. for each eventID, see if *all* columns match
        rows_identical = matches.all(axis=1)

        # 6. compute your counts
        total_in_both = len(common_ids)
        totally_unchanged = rows_identical.sum()
        changed_some_column = total_in_both - totally_unchanged
        percentage_of_changed_columns = changed_some_column / total_in_both

        return changed_some_column, total_in_both, totally_unchanged, percentage_of_changed_columns

    def sample_IM_metrics(self, path, sample_ratio, sampling_algo="AB", connectivity_threshold=0.5):

        # import log
        Evaluation.import_ocel(self,path)

        # my sampling
        sampling = smpl.OCEL_sampling()
        self.sample_pm4py = sampling.striped_sampling(self.ocel_pm4py, self.ocel_ocpa, sample_ratio, sampling_algo, connectivity_threshold)
        Evaluation.update_sample_from_pm4py_to_ocpa(self)

        #IM
        # ocel = ocel_import_factory.apply('result.xmlocel')
        print("reults-xml imported")
        ocpn = ocpn_discovery_factory.apply(self.sample_ocpa, parameters={"debug": True})
        print("Petri net mined")
        # metrics on model and sample

        Evaluation.calculate_quality_metrics(self, self.sample_ocpa, ocpn)

        print(self.quality_metrics)

    def eval_from_file(self, path_log_for_model, path_log_for_eval, token):

        #IM
        self.token_replay = token
        self.sample_ocpa = ocel_import_factory.apply(path_log_for_model)
        print("reults-xml imported")
        ocpn = ocpn_discovery_factory.apply(self.sample_ocpa, parameters={"debug": True})
        os.environ["PATH"] += os.pathsep + 'C:/Program Files/Graphviz/bin/'
        ocpn_vis_factory.save(ocpn_vis_factory.apply(ocpn), "oc_petri_net.png")
        print("Petri net mined")
        # metrics on model and sample

        self.ocel_ocpa_orginal = ocel_import_factory.apply(path_log_for_eval)
        start_time = time.time()
        Evaluation.calculate_quality_metrics(self, self.ocel_ocpa_orginal, ocpn)
        end_time = time.time()
        print("Time for eval was: " + str((end_time - start_time)) + "seconds ")
        print(self.quality_metrics)

        data = {'Info': ['path_log_for_model','path_log_for_eval','token','metrics' ],
                'Value': [path_log_for_model, path_log_for_eval, token, self.quality_metrics]}
        df = pd.DataFrame(data)
        df.to_csv("output.csv", index=False)


