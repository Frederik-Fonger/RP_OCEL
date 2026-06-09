import sampling.manager
import evaluation.OCEL_evaluation



if __name__ == '__main__':

    # -------------- Sampling -----------------

    sampling = sampling.manager.OCEL_sampling()
    sampling.apply('data/lrms collection/02_p2p.xml', 0.1, sampling_algo="RP", connectivity_threshold=0.5, file_type="XML", create_folder_with_sample_and_meta_data=True)
    

    # --------------  Evaluation  --------------

    # eval = evaluation.OCEL_evaluation.Evaluation()
    
    # ------ fitness -----
    # eval.eval_from_file(path_log_for_model='data/lrms collection/02_p2p.xml',path_log_for_eval= 'data/lrms collection/RP_mod_sampling/RP_0_4/02_p2p_sample_RP_0.4.xmlocel', token=False)
    

    # ------ MAE and coverage -----
    # eval.calculate_MAE_and_coverage(sample="data/lrms collection/RP_mod_sampling/RP_0_5/02_p2p_sample_RP_0.5.xmlocel", original_log="data/lrms collection/02_p2p.xml")
