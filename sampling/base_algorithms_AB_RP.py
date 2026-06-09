from collections import defaultdict, Counter
import json
from copy import deepcopy
from typing import Dict

import numpy as np
import pm4py


def ocel_allbehaviour_sampling(Sampling):
    """

    This funkction applies AllBehavior sampling on the OCEL. Therefore, the sample ratio is used. The goal is the first step is to represend all DFRs at leat one in the sample. In the second step the sample is filled with Remainder Plus sampling until the desired sample size is reached.

    :param Sampling: The Smample-Object from the manager.py, that contains information over the OCEL
    :return: a list of all Process Enecuation that are selected for the sample.
    """


    # prepare variables like it is in original AB sampling

    # Todo gleich PEs zusammenfassen und häufigkeit festhalten
    pe_string_list, pe_list = prepare_PEs_for_RP(Sampling)

    # calculate the expected count of event for the sample
    Sampling.target_event_count_in_sample = int(Sampling.event_count_origial_ocel * Sampling.sampling_ratio)

    # sort variants based on their frequency in descending order
    pe_list_sorted = sorted(
        pe_string_list, key=lambda element: element[1]['count'], reverse=True)

    num_of_traces_in_ol = len(Sampling.pe_list_DFRs)

    # add a expected occurrence count for each pe
    for variant in pe_list_sorted:
        variant[1]['expected_occurrence'] = variant[1]['count'] * Sampling.sampling_ratio

    # Take a pe from the variant with the highest frequency (will late be used for initializing the sample)
    pe_num_of_most_commen_pe = deepcopy(pe_list_sorted[0][1]['indices'][0])

    # delete the before added pe from the pe list and decrease the count of this variant by 1
    pe_list_sorted[0][1]['indices'].pop(0)
    pe_list_sorted[0][1]['count'] -= 1
    pe_list_sorted[0][1]['expected_occurrence'] -= 1

    # if the before added pe was the only one of the variant, then remove the pe-variant from the list
    if pe_list_sorted[0][1]['count'] == 0:
        pe_list_sorted.pop(0)

    # keep track of the variants that have been already sampled
    variants_that_have_been_sampled = []
    variants_that_have_been_sampled.append(
        pe_list[pe_num_of_most_commen_pe])

    # add the first variant to the sample with frequency 1
    allbehaviour_sample = []
    allbehaviour_sample.append(pe_num_of_most_commen_pe)

    # Determine the unsampled behaviour...Therefore check which behaviour is not part of our sample yet.
    metrics = _calculate_ratios(
        Sampling, sampled_event_log=allbehaviour_sample, pe_list=pe_list)

    # The unsampled behaviour list is a list of behaviour pairs that are part of the event log and haven't been added to the sample yet.
    unsampled_behaviour = [tuple(x[0])
                           for x in metrics.get("unsampled_behavior_list")]



    while calculate_size_of_sample(Sampling, allbehaviour_sample) < Sampling.target_event_count_in_sample and len(unsampled_behaviour) > 0:

        index_of_new_variant = 0
        max_normalized_count = 0

        # Iterate through every variant and determine the variant's dfr.
        # For each variant count the number of dfr that haven't been added to the sample yet.
        # Due to the fact that long variants (number of events in the sequence/trace is high) have a higher probability of having a high count (a high number of unsampled behaviour pairs),
        # divide the count of unsampled dfr in the variant by the count of dfr in the variant (normalization)
        for idx, variant in enumerate(pe_list_sorted):

            # get the object centric DFRs in the PE (variant)
            if(variant[1]['count']==0):
                print(variant)
            pairs_in_one_variant = pe_list[variant[1]['indices'][0]]
            pairs_in_one_variant = [x[0] for x in pairs_in_one_variant]

            # todo test if this counter works
            count_of_behaviour_pairs_that_are_part_of_the_variant_but_not_sampled_yet = 0

            # count how many of the DFR in this pe/variant are unsampled
            for pair in pairs_in_one_variant:
                if pair in unsampled_behaviour:
                    count_of_behaviour_pairs_that_are_part_of_the_variant_but_not_sampled_yet += 1

            # normalization the number of unsampled DFRs in this pe/variant
            if count_of_behaviour_pairs_that_are_part_of_the_variant_but_not_sampled_yet > 0:
                count_of_behaviour_pairs_that_are_part_of_the_variant_but_not_sampled_yet /= len(
                    pairs_in_one_variant)
            else:
                count_of_behaviour_pairs_that_are_part_of_the_variant_but_not_sampled_yet = -1

            # We need to safe the index of the variant with the highest normalized count.
            if count_of_behaviour_pairs_that_are_part_of_the_variant_but_not_sampled_yet > max_normalized_count:
                max_normalized_count = count_of_behaviour_pairs_that_are_part_of_the_variant_but_not_sampled_yet
                index_of_new_variant = idx

        # get the index of one pe from the variant that is selected to be added to the sample
        pe_num_to_add_to_sample = pe_list_sorted[index_of_new_variant][1]['indices'][0]
        # remove th index of the pe from the list of pe´s from the variant
        pe_list_sorted[index_of_new_variant][1]['indices'].pop(0)
        # As one pe was removed in the line above,  the pe count in the variant need to be decreases by one.
        pe_list_sorted[index_of_new_variant][1]['count'] -= 1
        # As the is now one pe of this variant added to the sample, the expected (remaining) occurrence is also decreased by one.
        pe_list_sorted[index_of_new_variant][1]['expected_occurrence'] -= 1
        # if the before added pe was the only one of the variant, then remove the pe-variant from the list
        if pe_list_sorted[index_of_new_variant][1]['count'] == 0:
            pe_list_sorted.pop(index_of_new_variant)

        # The variant with the highest normalized count will be added to our sample and...
        allbehaviour_sample.append(pe_num_to_add_to_sample)
        # Update the list of already sampled DRFs
        variants_that_have_been_sampled.append(pe_list[pe_num_to_add_to_sample])

        # we remove the behaviour of the new variant (variant with the highest normalized count) from our unsampled behaviour list
        #                                               v HERE: the pe_list is a list of tupels with the DFRs and the count of the pe. Here needed are only the DFRs.
        behaviour_of_the_new_variant_to_be_removed = [x[0] for x in pe_list[pe_num_to_add_to_sample]]

        # update the list of unsampled DFRs
        unsampled_behaviour = list((set(unsampled_behaviour)).difference(
            set(behaviour_of_the_new_variant_to_be_removed)))

        # print("Size limit: " + str(calculate_size_of_sample(Sampling, allbehaviour_sample)) + " < " + str(Sampling.target_event_count_in_sample) + " and Lenght us. B: " + str(len(unsampled_behaviour)))
        # print("test")
    # We have all behaviour pairs of the event log in our sample. In case we have remaining free slots in our sample we try to improve the sample's representativeness by using the remainderplus algorithm

    # update the sorting of the pe_list
    pe_list_sorted = sorted(
        pe_list_sorted, key=lambda element: element[1]['expected_occurrence'], reverse=True)

    # Check for every variant if the expected occurrenc is over 1 and add the rounded down amount of pe´s of this variant to the sample.
    for variant in pe_list_sorted:

        # append the variant to the sample if the variants expected occurrence is greater or equal to 1.
        # The sampled variant's occurrence is reassigned with the integer part of the calculated variants expected occurrence.
        if variant[1]['expected_occurrence'] >= 1:

            # calculated th rounded down amount of expected occurrences of this variant
            intnum = int(variant[1]['expected_occurrence'])
            # calcualted the deciamt places that are not used by the aforementioned rounding down.
            remainder = variant[1]['expected_occurrence'] % intnum
            variant[1]['expected_occurrence'] = remainder

            # create a list of pe´s from this variant
            corpus = []
            corpus = corpus + (variant[1]['indices'][:intnum])
            allbehaviour_sample.extend(corpus)


    #################From here same as RemainderPlus-Sampling#############################

 # we have to sort on behaviour characteristics in case two or more PEs have the same remainder. Therefore we first check which behaviour is undersampled or oversampled
    intermediate_results = _calculate_ratios(Sampling, sampled_event_log=allbehaviour_sample, pe_list=pe_list)

    behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version = {tuple(
        x): y for x, y in intermediate_results.get("list_of_pairs_with_sample_ratio")}
    behaviour_pairs_in_original_log_with_count = {
        tuple(x): y for x, y in intermediate_results.get("list_of_pairs_with_count_ol")}

    # initiliaze every variant with rank = 0
    for variant in pe_list_sorted:
        variant.append(0)


    while calculate_size_of_sample(Sampling, allbehaviour_sample) < Sampling.target_event_count_in_sample:

        # Calculate each PE's rank (it's called "differenzwert")
        # The PE will get a positive rank if the variant has more undersampled behaviour than oversampled behaviour and will likely be sampled
        # The rank of the PE will increase by one for each undersampled behaviour and will decrease by one for each oversampled behaviour in the sample
        for variant in pe_list_sorted:

            rank = 0

            # get the object centric DFRs in the PE
            pairs_in_one_variant = pe_list[variant[1]['indices'][0]]
            pairs_in_one_variant = [x[0] for x in pairs_in_one_variant]

            for pair in pairs_in_one_variant:

                if behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version.get(pair) > Sampling.sampling_ratio:
                    # pair is oversampled
                    rank -= 1
                if behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version.get(pair) < Sampling.sampling_ratio:
                    # pair is undersampled
                    rank += 1
                else:
                    # pair is perfectly sampled
                    rank += 0

            # normalize rank by the number of pairs in the variant
            if (len(pairs_in_one_variant) > 0):
                rank = rank / len(pairs_in_one_variant)

            # assign the rank to the variant
            variant[2] = rank

        # sort by two attributes. First by the remainder and then sort by the rank
        pe_list_sorted = sorted(
            pe_list_sorted, key=lambda x: (x[1]['expected_occurrence'], x[2]), reverse=True)

        # the variant that has to sampled next has to be on index 0 of the list (highest remainder and compared to all variants with the same remainder it has the highest rank)
        variant_to_be_sampled = pe_list_sorted[0]

        # update all corrsponding behaviour pairs:
        variants_pairs_to_be_sampled = pe_list[variant_to_be_sampled[1]['indices'][0]]
        variants_pairs_to_be_sampled = [x[0] for x in variants_pairs_to_be_sampled]

        for pair in variants_pairs_to_be_sampled:
            behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version[pair] = ((
                                                                                               behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version.get(
                                                                                                   pair) * behaviour_pairs_in_original_log_with_count.get(
                                                                                           pair)) + 1) / behaviour_pairs_in_original_log_with_count.get(
                pair)

        # Add the variant to the sample. Therefore we have to check whether the variant is already part of our sample. In case the variant is already in our sample, we have to update the frequency of that variant in our sample,
        # Otherwise we add the variant to our sample with the frequency of one
        corpus = []
        corpus.append(variant_to_be_sampled[1]['indices'][0])
        allbehaviour_sample.extend(corpus)


        # update the remainder because we have sampled the unique variant now
        variant_to_be_sampled[1]['expected_occurrence'] = 0

    # return selected PEs
    return allbehaviour_sample



def ocel_remainder_plus_sampling(Sampling):

    # prepare variables like it is in original RP sampling

    # Todo gleich PEs zusammenfassen und häufigkeit festhalten
    pe_string_list, pe_list = prepare_PEs_for_RP(Sampling)

    # calcuate the tages count of event for the sample
    Sampling.target_event_count_in_sample = int(Sampling.event_count_origial_ocel * Sampling.sampling_ratio)

    # sort variants based on their frequency in descending order
    pe_list_sorted = sorted(
        pe_string_list, key=lambda element: element[1]['count'], reverse=True)



    # In the beginning no trace is part of the sample. We build up the sample successively
    remainder_plus_sample = []

    # iterate through every PE in the event log...
    for variant in pe_list_sorted:
        # and calculate the PE' expected occurrences in the sample
        the_variants_expected_occurrence = variant[1]['count'] * Sampling.sampling_ratio
        variant[1]['expected_occurrence'] = the_variants_expected_occurrence

        # append the PE to the sample if the variants expected occurrence is greater or equal to 1.
        # The sampled PE's occurrence is reassigned with the integer part of the calculated variants expected occurrence.
        if variant[1]['expected_occurrence'] >= 1:
            intnum = int(variant[1]['expected_occurrence'])
            remainder = variant[1]['expected_occurrence'] % intnum
            variant[1]['expected_occurrence'] = remainder
            corpus = []
            corpus = corpus + (variant[1]['indices'][:intnum])
            remainder_plus_sample.extend(corpus)


    # we have to sort on behaviour characteristics in case two or more PEs have the same remainder. Therefore we first check which behaviour is undersampled or oversampled
    intermediate_results = _calculate_ratios(Sampling, sampled_event_log=remainder_plus_sample, pe_list=pe_list)

    behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version = {tuple(
        x): y for x, y in intermediate_results.get("list_of_pairs_with_sample_ratio")}
    behaviour_pairs_in_original_log_with_count = {
        tuple(x): y for x, y in intermediate_results.get("list_of_pairs_with_count_ol")}

    # initiliaze every variant with rank = 0
    for variant in pe_list_sorted:
        variant.append(0)


    while calculate_size_of_sample(Sampling, remainder_plus_sample) < Sampling.target_event_count_in_sample:

        # Calculate each PE's rank (it's called "differenzwert")
        # The PE will get a positive rank if the variant has more undersampled behaviour than oversampled behaviour and will likely be sampled
        # The rank of the PE will increase by one for each undersampled behaviour and will decrease by one for each oversampled behaviour in the sample
        for variant in pe_list_sorted:

            rank = 0

            # get the object centric DFRs in the PE
            pairs_in_one_variant = pe_list[variant[1]['indices'][0]]
            pairs_in_one_variant = [x[0] for x in pairs_in_one_variant]

            for pair in pairs_in_one_variant:

                if behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version.get(pair) > Sampling.sampling_ratio:
                    # pair is oversampled
                    rank -= 1
                if behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version.get(pair) < Sampling.sampling_ratio:
                    # pair is undersampled
                    rank += 1
                else:
                    # pair is perfectly sampled
                    rank += 0

            # normalize rank by the number of pairs in the variant
            if (len(pairs_in_one_variant) > 0):
                rank = rank / len(pairs_in_one_variant)

            # assign the rank to the variant
            variant[2] = rank

        # sort by two attributes. First by the remainder and then sort by the rank
        pe_list_sorted = sorted(
            pe_list_sorted, key=lambda x: (x[1]['expected_occurrence'], x[2]), reverse=True)

        # the variant that has to sampled next has to be on index 0 of the list (highest remainder and compared to all variants with the same remainder it has the highest rank)
        variant_to_be_sampled = pe_list_sorted[0]

        # update all corrsponding behaviour pairs:
        variants_pairs_to_be_sampled = pe_list[variant_to_be_sampled[1]['indices'][0]]
        variants_pairs_to_be_sampled = [x[0] for x in variants_pairs_to_be_sampled]

        for pair in variants_pairs_to_be_sampled:
            behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version[pair] = ((
                                                                                               behaviour_pairs_with_behaviour_pair_sample_ratio_dict_version.get(
                                                                                                   pair) * behaviour_pairs_in_original_log_with_count.get(
                                                                                           pair)) + 1) / behaviour_pairs_in_original_log_with_count.get(
                pair)

        # Add the variant to the sample. Therefore we have to check whether the variant is already part of our sample. In case the variant is already in our sample, we have to update the frequency of that variant in our sample,
        # Otherwise we add the variant to our sample with the frequency of one
        corpus = []
        corpus.append(variant_to_be_sampled[1]['indices'][0])
        remainder_plus_sample.extend(corpus)


        # update the remainder because we have sampled the unique variant now
        variant_to_be_sampled[1]['expected_occurrence'] = 0

    # return selected PEs
    return remainder_plus_sample



def prepare_PEs_for_RP(Sampling):
    """
    This function prepares the format of the process execution for the following sampling.

    :param Sampling: OCEL Sampling object
    :return:    pe_string_list: list of PEs with pe as stings
                pe_list:        list of pe with indexes
    """

    strings = []
    pe_list = []


    for pe in Sampling.pe_list_DFRs:

        # create a list of the PEs with a string of the activties
        strings.append(str(pe.items()))
        # create a list of the PEs with the DFRs from the PEs and the corresponding count of the DFR
        pe_list.append(list(pe.items()))

        # Dictionary to store counts and indices
    count_indices = defaultdict(lambda: {"count": 0, "indices": []})

    # Process the list
    for index, s in enumerate(strings):
        count_indices[s]["count"] += 1
        count_indices[s]["indices"].append(index)

    # Convert defaultdict to a regular dictionary
    count_of_same_PEs = dict(count_indices)
    pe_string_list = []
    for key, value in count_of_same_PEs.items():
        temp = [key]
        temp.append(value)
        pe_string_list.append(temp)

    return pe_string_list, pe_list


def _calculate_ratios(Sampling, sampled_event_log, pe_list):

    """
       This function calculates the ratios of behavior pairs between an original event log and a sampled event log.
       It creates a ratio matrix to compare the frequencies of behavior pairs in both logs.
       This function is used for intermediate calculations and should not be used directly by the user.

    Parameters:


    Returns:
        unsampled_behavior_list: A list of behavior pairs that were not sampled.
        list_of_pairs_with_sample_ratio: A list of behavior pairs along with their sampling ratios.
        list_of_pairs_with_count_ol: The frequency of different behavior pairs in the original log.
    """

    # count the frequency of the different pairs for the original log
    list_of_pairs_with_count_ol = flatten_pe_list(pe_list)

    pe_in_sample = [pe_list[i] for i in sampled_event_log]

    # count the frequency of the different pairs for the sample log
    list_of_pairs_with_count_sl = flatten_pe_list(pe_in_sample)

    # extract list of unique activities
    list_of_single_unique_activities = _get_list_of_single_unique_activities(
        list_of_pairs_with_count_ol)

    # Sort activities
    list_of_single_unique_activities_sorted = sorted(
        [x for x in list_of_single_unique_activities])

    # initialize adjacency matrix:
    # number of rows and columns = number of unique activities in the original event log
    num_of_rows_and_columns = len(list_of_single_unique_activities)

    adjacency_matrix = np.zeros(
        (num_of_rows_and_columns, num_of_rows_and_columns), dtype=float)

    for behavior_with_count in list_of_pairs_with_count_ol:
        behavior_start = behavior_with_count[0][0]
        behavior_end = behavior_with_count[0][1]
        behavior_count = behavior_with_count[1]

        for idx, element in enumerate(list_of_single_unique_activities_sorted):
            if element == behavior_start:
                row = idx
            if element == behavior_end:
                col = idx

        adjacency_matrix[row][col] = behavior_count

    adjacency_matrix_sample = np.zeros(
        (num_of_rows_and_columns, num_of_rows_and_columns), dtype=float)

    for behavior_with_count in list_of_pairs_with_count_sl:
        behavior_start = behavior_with_count[0][0]
        behavior_end = behavior_with_count[0][1]
        behavior_count = behavior_with_count[1]

        for idx, element in enumerate(list_of_single_unique_activities_sorted):
            if element == behavior_start:
                row = idx
            if element == behavior_end:
                col = idx

        adjacency_matrix_sample[row][col] = behavior_count

    ################## build DFR-Ratio Matrix#############################################
    adjacency_matrix_sample_ratio = np.zeros(
        (num_of_rows_and_columns, num_of_rows_and_columns), dtype=float)

    list_of_pairs_with_sample_ratio = [
        [list(x), y] for x, y in list_of_pairs_with_count_ol]

    for behavior_with_count in list_of_pairs_with_sample_ratio:
        behavior_start = behavior_with_count[0][0]
        behavior_end = behavior_with_count[0][1]
        behavior_count = behavior_with_count[1]

        for idx, element in enumerate(list_of_single_unique_activities_sorted):
            # identify the right column and row
            if element == behavior_start:
                row = idx

            if element == behavior_end:
                col = idx

        # caluculate the ratio of the dfr ...
        if adjacency_matrix_sample[row][col] <= behavior_count and adjacency_matrix_sample[row][col] > 0:
            # adjacency_matrix_sample_ratio[row][col] = round(adjacency_matrix_sample[row][col]/behavior_count,2)
            adjacency_matrix_sample_ratio[row][col] = adjacency_matrix_sample[row][col] / behavior_count
            behavior_with_count[1] = adjacency_matrix_sample_ratio[row][col]
        elif adjacency_matrix_sample[row][col] == 0:
            behavior_with_count[1] = adjacency_matrix_sample_ratio[row][col]
        else:
            pass

        unsampled_behavior_list = []

        for x in list_of_pairs_with_sample_ratio:

            if x[1] == 0:
                unsampled_behavior_list.append(x)

    calculated_ratios = {
        "unsampled_behavior_list": unsampled_behavior_list,
        "list_of_pairs_with_sample_ratio": list_of_pairs_with_sample_ratio,
        "list_of_pairs_with_count_ol": list_of_pairs_with_count_ol,

    }

    return calculated_ratios


def flatten_pe_list(pe_list):
    """
    Sums the counts of DFRs in the pe list. The calculate the total count of a DFR over all PEs. The count from all PEs
    of the DFRs needs to be sumed. Thes function dose this for all the DFRs in all the PEs in the pe_list.

    :param pe_list:
    :return:
    """
    list_of_pairs_with_count_ol = []
    for pe in pe_list:
        list_of_pairs_with_count_ol.extend(pe)

    # Create a default dict to store sums for each index tuple
    index_sums = defaultdict(int)

    # Loop through the list and sum numbers by index tuple
    for index, number in list_of_pairs_with_count_ol:
        index_sums[index] += number

    list_of_pairs_with_count_ol = []
    for key, value in index_sums.items():
        temp = [key]
        temp.append(value)
        list_of_pairs_with_count_ol.append(temp)

    return  list_of_pairs_with_count_ol

def _get_list_of_single_unique_activities(list_of_pairs):
    unique_activities_set = set()
    for a_tuple in list_of_pairs:
        unique_activities_set.update(a_tuple[0])
    return list(unique_activities_set)

def calculate_size_of_sample(Sampling, remainder_plus_sample):

    sample_size = 0
    for pe in remainder_plus_sample:
        sample_size += len(Sampling.ocel_ocpa.process_executions[pe])
    return sample_size


