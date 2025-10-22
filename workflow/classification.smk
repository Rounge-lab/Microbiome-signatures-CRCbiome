import csv
import pandas as pd

configfile: "config/ml_config.yaml"

model_list=pd.read_csv(config["model_list"], sep = '\t')

ml_algs=["lasso", "rf", "xgb", "svm", "nnet"]

feature_strat =  model_list[ model_list["model_type"] == "feature-strat"]["dataset"].to_list()

outcome_strat =  model_list[ model_list["model_type"] == "outcome-strat"]["dataset"].to_list()

full_datasets =  model_list[ model_list["model_type"] == "base"]["dataset"].to_list()

datasets =  feature_strat + outcome_strat + full_datasets

n_iter = config["n_iter"]
iters = [ str(i) for i in range(n_iter)]

min_feature_frac, max_feature_frac = config["min_feature_frac"], config["max_feature_frac"]
n_f = [n for n in [ *range(min_feature_frac,max_feature_frac)] if n % 2 == 0 ]

n_sample_fractions = config["n_sample_fractions"]
min_sample_fraction = config["min_sample_fraction"]
fractions = [ i+1 for i in range(n_sample_fractions)]
n_sample_reps = [i + 1 for i in range(5)]

n_features_datasets = full_datasets
n_samples_datasets = full_datasets

analysis_name="test"

final_datasets = ["mags", "mphl"]
secondary_lr_analysis_name = analysis_name

rule run_primary_models:
    input:
        "data/model_training/" + analysis_name + "/gathered_primary_probs.tsv",

rule run_n_features:
    input:
        "data/model_training/" + analysis_name + "/n_features_gathered_probs.tsv",

rule run_n_samples:
    input:
        "data/model_training/" + analysis_name + "/sample_frac_gathered_model_specs.tsv",

rule run_secondary_models:
    input:
        "data/model_training/" + analysis_name + "/gathered_secondary_assessment_probs.tsv",

rule run_final_primary_model:
    input:
        expand("data/models/final_model/{ml_alg}_{dataset}/probs.tsv", ml_alg = ["rf"], dataset = final_datasets),

rule run_final_models:
    input:
        expand("data/models/final_model/{ml_alg}_{dataset}/probs.tsv", ml_alg = ["rf"], dataset = final_datasets),
        ## Secondary models
        expand("data/models/final_model/{ml_alg}_{dataset}/secondary_lr/{preselection}/{add_var}_probs.tsv", 
                ml_alg = ["rf"], dataset = final_datasets,
                preselection=["norm_nn_v_aa_crc"], 
                add_var = ["fit", "fit_demo_wcrf"])

rule secondary_models:
    input:
        expand("data/model_training/" + analysis_name + "/secondary_lr/{preselection}/{m_type}-{add_var}_probs.tsv", 
                preselection = ["norm_nn_v_aa_crc"], 
                # m_type = ["base", "n_features", "sample_frac"],
                m_type = ["base"],
                add_var = ["fit", "fit_demo", "fit_demo_wcrf", "fit_lididem", "catfit", "catfit_demo", "catfit_demo_wcrf", "catfit_lididem", "restr_PPI", "full_PPI"]),
    output:
        "data/model_training/" + analysis_name + "/secondary_training.flag"
    shell:
        "touch {output}"

rule train_val_no_tune:
    input:
        training=config["dataset_folder"]+"{dataset}/train.tsv",
        assessment=config["dataset_folder"]+"{dataset}/assessment_set.tsv"
    output:
        model_specs="data/models/{ml_alg}_{dataset}/training/model_specs.tsv",
        probs="data/models/{ml_alg}_{dataset}/training/probs.tsv",
        assessment_probs="data/models/{ml_alg}_{dataset}/training/assessment_probs.tsv",
        feature_importance="data/models/{ml_alg}_{dataset}/training/feature_importance.tsv",
    params:
        iterations=iters,
        alg_params=lambda w: config["ml_params"][w.ml_alg],
        dataset_proc_steps=config["dataset_proc_steps"],
        target_var=config["target_var"],
        training_prop=0.8,
        bootstraps=50,
        folds=5,
        repeats=20,
        resampling="vfold_cv",
        var_info=config["dataset_folder"]+"{dataset}/var_info_dontuse.tsv"
    conda:
        "envs/tidymodels.yaml"
    threads:
        10
    script:
        "scripts/classification/ml_resampling_wrapper.R"

rule sec_mod:
    input:
        training="data/models/{prim_alg}_{dataset}/training/probs.tsv",
        assessment="data/models/{prim_alg}_{dataset}/training/assessment_probs.tsv",
        add_vars_train=config["dataset_folder"]+"secondary_lr/{preselection}/{add_var}/train.tsv",
        add_vars_assess=config["dataset_folder"]+"secondary_lr/{preselection}/{add_var}/assessment_set.tsv",
    output:
        model_specs="data/models/{prim_alg}_{dataset}/training_sec/{preselection}/{add_var}/{sec_alg}_model_specs.tsv",
        probs="data/models/{prim_alg}_{dataset}/training_sec/{preselection}/{add_var}/{sec_alg}_probs.tsv",
        assessment_probs="data/models/{prim_alg}_{dataset}/training_sec/{preselection}/{add_var}/{sec_alg}_assessment_probs.tsv",
        feature_importance="data/models/{prim_alg}_{dataset}/training_sec/{preselection}/{add_var}/{sec_alg}_feature_importance.tsv",
    log:
        log="logs/modelling/{prim_alg}_{dataset}/{preselection}/{add_var}/{sec_alg}_log.txt"
    params:
        alg_params=lambda w: config["ml_params"][w.sec_alg],
        dataset_proc_steps=config["dataset_proc_steps"],
        target_var=config["target_var"],
        resampling="predefined",
        pred=".pred_positive"
    conda:
        "envs/tidymodels.yaml"
    threads:
        1
    script:
        "scripts/classification/ml_resampling_wrapper.R"

def sec_pred_output_files(w, mod_file):
    files = expand("data/models/"+w.prim_alg+"_"+w.dataset+"/training_sec/{preselection}/{add_var}/{sec_alg}_"+mod_file+".tsv",
        preselection = ["norm_nn_v_aa_crc"],
        add_var = ["fit", "fit_demo", "fit_demo_wcrf", "fit_lididem", 
                    "catfit", "catfit_demo", "catfit_demo_wcrf", "catfit_lididem", 
                    "restr_PPI", "full_PPI"],
        sec_alg = ["rf", "xgb", "svm", "nnet", "lasso", "log_reg"])
    return files

localrules: gather_sec_mod_res

rule gather_sec_mod_res:
    input:
        model_specs = lambda w: sec_pred_output_files(w, "model_specs"),
        probs = lambda w: sec_pred_output_files(w, "probs"),
        assessment_probs = lambda w: sec_pred_output_files(w, "assessment_probs"),
        feature_importance = lambda w: sec_pred_output_files(w, "feature_importance")
    output:
        gathered_model_specs="data/models/{prim_alg}_{dataset}/training_sec/gathered_model_specs.tsv",
        gathered_assessment="data/models/{prim_alg}_{dataset}/training_sec/gathered_assessment_probs.tsv",
        gathered_probs="data/models/{prim_alg}_{dataset}/training_sec/gathered_probs.tsv",
        gathered_fi="data/models/{prim_alg}_{dataset}/training_sec/gathered_fi.tsv"
    params:
        prim_mod = "{prim_alg}_{dataset}",
        missing_col = "no_m"
    threads:
        1
    script:
        "scripts/classification/gather_no_tune_training.py"

def input_for_gather(mod_file, stage):
    if stage == "primary":
        if mod_file == "fi":
            mod_file = "feature_importance"
        return [ "data/models/" + row["alg"] + "_" + row["dataset"] + "/training/"+mod_file+".tsv" 
                                for _, row in model_list[ model_list["dataset"].isin(datasets)].iterrows()]
    if stage == "secondary":
        return [ "data/models/" + row["alg"] + "_" + row["dataset"] + "/training_sec/gathered_"+mod_file+".tsv" 
                                for _, row in model_list[ model_list["dataset"].isin(datasets)].iterrows()
                                if not "neg-v-serr" in row["dataset"]]

localrules: gather_training_res

rule gather_training_res:
    input:
        model_specs = lambda w: input_for_gather("model_specs", w.stage),
        probs = lambda w: input_for_gather("probs", w.stage),
        assessment_probs = lambda w: input_for_gather("assessment_probs", w.stage),
        feature_importance = lambda w: input_for_gather("fi", w.stage),
    output:
        gathered_model_specs="data/model_training/" + analysis_name + "/gathered_{stage}_model_specs.tsv",
        gathered_assessment="data/model_training/" + analysis_name + "/gathered_{stage}_assessment_probs.tsv",
        gathered_probs="data/model_training/" + analysis_name + "/gathered_{stage}_probs.tsv",
        gathered_fi="data/model_training/" + analysis_name + "/gathered_{stage}_fi.tsv"
    threads:
        1
    script:
        "scripts/classification/gather_no_tune_training.py"

rule n_features:
    input:
        training=config["dataset_folder"]+"{dataset}/train.tsv",
        assessment=config["dataset_folder"]+"{dataset}/assessment_set.tsv"
    output:
        model_specs=temp("data/models/{ml_alg}_{dataset}/training/model_specs_{n_features}.tsv"),
        probs=temp("data/models/{ml_alg}_{dataset}/training/probs_{n_features}.tsv"),
        assessment_probs=temp("data/models/{ml_alg}_{dataset}/training/assessment_probs_{n_features}.tsv"),
        feature_importance=temp("data/models/{ml_alg}_{dataset}/training/feature_importance_{n_features}.tsv"),
    params:
        iterations=iters,
        alg_params=lambda w: config["ml_params"][w.ml_alg],
        dataset_proc_steps=config["dataset_proc_steps"],
        target_var=config["target_var"],
        training_prop=0.8,
        bootstraps=50,
        folds=5,
        repeats=20,
        resampling="vfold_cv",
        var_info=config["dataset_folder"]+"{dataset}/var_info.tsv",
    conda:
        "envs/tidymodels.yaml"
    threads:
        10
    script:
        "scripts/classification/ml_resampling_wrapper.R"

localrules: gather_n_features

rule gather_n_features:
    input:
        model_specs=expand("data/models/{ml_alg}_{dataset}/training/model_specs_{n_features}.tsv", 
                                ml_alg = ml_algs, dataset=n_features_datasets, n_features = n_f),
        probs=expand("data/models/{ml_alg}_{dataset}/training/probs_{n_features}.tsv", 
                                ml_alg = ml_algs, dataset=n_features_datasets, n_features = n_f),
        assessment_probs=expand("data/models/{ml_alg}_{dataset}/training/assessment_probs_{n_features}.tsv", 
                                ml_alg = ml_algs, dataset=n_features_datasets, n_features = n_f),
        feature_importance=expand("data/models/{ml_alg}_{dataset}/training/feature_importance_{n_features}.tsv", 
                                ml_alg = ml_algs, dataset=n_features_datasets, n_features = n_f)
    output:
        gathered_model_specs="data/model_training/" + analysis_name + "/n_features_gathered_model_specs.tsv",
        gathered_assessment="data/model_training/" + analysis_name + "/n_features_gathered_assessment_probs.tsv",
        gathered_probs="data/model_training/" + analysis_name + "/n_features_gathered_probs.tsv",
        gathered_fi="data/model_training/" + analysis_name + "/n_features_gathered_fi.tsv"
    threads:
        1
    script:
        "scripts/classification/gather_no_tune_training.py"

rule sample_frac:
    input:
        training=config["dataset_folder"]+"{dataset}/train.tsv",
        assessment=config["dataset_folder"]+"{dataset}/assessment_set.tsv"
    output:
        model_specs=temp("data/models/{ml_alg}_{dataset}/training_sample_frac/model_specs_{fraction}_{rep}.tsv"),
        probs=temp("data/models/{ml_alg}_{dataset}/training_sample_frac/probs_{fraction}_{rep}.tsv"),
        assessment_probs=temp("data/models/{ml_alg}_{dataset}/training_sample_frac/assessment_probs_{fraction}_{rep}.tsv"),
        feature_importance=temp("data/models/{ml_alg}_{dataset}/training_sample_frac/feature_importance_{fraction}_{rep}.tsv"),
    params:
        iterations=iters,
        alg_params=lambda w: config["ml_params"][w.ml_alg],
        dataset_proc_steps=config["dataset_proc_steps"],
        target_var=config["target_var"],
        folds=5,
        repeats=1,
        resampling="vfold_cv",
        var_info=config["dataset_folder"]+"{dataset}/var_info.tsv",
        min_sample_frac=min_sample_fraction,
        sample_fracs=len(fractions)
    conda:
        "envs/tidymodels.yaml"
    threads:
        5
    script:
        "scripts/classification/ml_resampling_wrapper.R"

localrules: gather_sample_frac

rule gather_sample_frac:
    input:
        model_specs=expand("data/models/{ml_alg}_{dataset}/training_sample_frac/model_specs_{fraction}_{rep}.tsv", 
                                ml_alg = ml_algs, dataset=n_samples_datasets, fraction = fractions, rep=n_sample_reps),
        probs=expand("data/models/{ml_alg}_{dataset}/training_sample_frac/probs_{fraction}_{rep}.tsv", 
                                ml_alg = ml_algs, dataset=n_samples_datasets, fraction = fractions, rep=n_sample_reps),
        assessment_probs=expand("data/models/{ml_alg}_{dataset}/training_sample_frac/assessment_probs_{fraction}_{rep}.tsv", 
                                ml_alg = ml_algs, dataset=n_samples_datasets, fraction = fractions, rep=n_sample_reps),
        feature_importance=expand("data/models/{ml_alg}_{dataset}/training_sample_frac/feature_importance_{fraction}_{rep}.tsv", 
                                ml_alg = ml_algs, dataset=n_samples_datasets, fraction = fractions, rep=n_sample_reps)
    output:
        gathered_model_specs="data/model_training/" + analysis_name + "/sample_frac_gathered_model_specs.tsv",
        gathered_assessment="data/model_training/" + analysis_name + "/sample_frac_gathered_assessment_probs.tsv",
        gathered_probs="data/model_training/" + analysis_name + "/sample_frac_gathered_probs.tsv",
        gathered_fi="data/model_training/" + analysis_name + "/sample_frac_gathered_fi.tsv"
    threads:
        1
    script:
        "scripts/classification/gather_no_tune_training.py"


def get_model_files(wildcards, filetype):
    
    pre={"base": "/", "sample_frac": "/sample_frac_", "n_features": "/n_features_"}
    
    if filetype == "probs":
        return "data/model_training/"+analysis_name+pre.get(wildcards.m_type)+"gathered_probs.tsv"
    if filetype == "assess":
        return "data/model_training/"+analysis_name+pre.get(wildcards.m_type)+"gathered_assessment_probs.tsv"


rule riskfac_secondary_logreg:
    input:
        probs=lambda wc: get_model_files(wc, "probs"),
        assessment_probs=lambda wc: get_model_files(wc, "assess"),
        add_vars_train=config["dataset_folder"]+"secondary_lr/{preselection}/{add_var}/train.tsv",
        add_vars_assess=config["dataset_folder"]+"secondary_lr/{preselection}/{add_var}/assessment_set.tsv",
    output:
        probs="data/model_training/" + analysis_name + "/secondary_lr/{preselection}/{m_type}-{add_var}_probs.tsv",
        assessment_probs="data/model_training/" + analysis_name + "/secondary_lr/{preselection}/{m_type}-{add_var}_assessment_probs.tsv",
        fi="data/model_training/" + analysis_name + "/secondary_lr/{preselection}/{m_type}-{add_var}_fi.tsv",
    params:
        steps=config["ml_params"]["secondary_lr"]["steps"],
        pred=".pred_positive",
        target="target",
        ref_level="negative"
    conda:
        "envs/tidymodels.yaml"
    threads:
        5
    script:
        "scripts/classification/ml_wrap_add_var_test_2.R"


rule final_model:
    input:
        train=config["dataset_folder"]+"{dataset}/train.tsv",
        test=config["dataset_folder"]+"{dataset}/test.tsv"
    output:
        model="data/models/final_model/{ml_alg}_{dataset}/model.Rds",
        model_specs="data/models/final_model/{ml_alg}_{dataset}/model_specs.tsv",
        probs="data/models/final_model/{ml_alg}_{dataset}/probs.tsv",
        feature_importance="data/models/final_model/{ml_alg}_{dataset}/feature_importance.tsv",
    params:
        alg_params=lambda w: config["ml_params"][w.ml_alg],
        dataset_proc_steps=config["dataset_proc_steps"],
        target_var=config["target_var"],
        resampling="none",
    conda:
        "envs/tidymodels.yaml"
    script:
        "scripts/classification/ml_test_set_wrapper.R"

rule final_secondary:
    input:
        train="data/models/final_model/{ml_alg}_{dataset}/probs.tsv",
        add_vars_train=config["dataset_folder"]+"secondary_lr/{preselection}/{add_var}/train.tsv",
        training_set_lr_parameters="data/model_training/"+secondary_lr_analysis_name+"/gathered_secondary_fi.tsv"
    output:
        secondary_predictions="data/models/final_model/{ml_alg}_{dataset}/secondary_lr/{preselection}/{add_var}_probs.tsv",
    params:
        steps=config["ml_params"]["secondary_lr"]["steps"],
        pred=".pred_positive",
        target="target",
        ref_level="negative"
    conda:
        "envs/tidymodels.yaml"
    threads:
        1
    script:
        "scripts/classification/secondary_models_final.R"


