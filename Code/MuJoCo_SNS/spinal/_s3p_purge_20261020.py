import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
for db, name in (("sqlite:///optuna_w2lvar.db", "curr_w2lvar_s3p_perturb"),
                 ("sqlite:///optuna_syn6.db", "curr_syn6_s3p_perturb")):
    try:
        optuna.delete_study(study_name=name, storage=db)
        print("deleted", name)
    except Exception as e:
        print("skip", name, e)
