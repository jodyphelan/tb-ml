import tb_ml
import pathlib


test_data = f"{pathlib.Path(__file__).resolve().parent}/test_data"
test_reads = f"{test_data}/test_aligned_reads.cram"


def run_test(capfd, cont_arg_list, res_file):
    # run the containers and capture the output
    tb_ml.run_containers(cont_arg_list)
    out, err = capfd.readouterr()
    # make sure there was no error msg
    assert err == ""
    # now test the output --> start with the header
    relevant_header = [line for line in out if line.startswith("# Container")]
    for line, (cont, args) in zip(relevant_header, cont_arg_list):
        arg_string = " ".join(args)
        assert line.strip() == f'# Container: {cont}; Args: "{arg_string}"'
    # now the actual results
    results = [line for line in out.strip().split("\n") if not line.startswith("#")]
    with open(res_file) as f:
        for line, res in zip(f.readlines(), results):
            assert line.strip() == res.strip()


def test_one_hot_into_neural_net(capfd):
    one_hot_container = "julibeg/tb-ml-one-hot-encoded-seqs-from-aligned-reads:v0.4.0"
    neural_net_container = (
        "julibeg/tb-ml-neural-net-from-one-hot-encoded-seqs-13-drugs:v0.7.0"
    )
    cont_arg_list = [
        [
            neural_net_container,
            "--get-target-loci -o target-loci.csv".split(),
        ],
        [
            one_hot_container,
            f"-b {test_reads} -r target-loci.csv -o one-hot-seqs.csv".split(),
        ],
        [
            neural_net_container,
            ["one-hot-seqs.csv"],
        ],
    ]
    run_test(capfd, cont_arg_list, f"{test_data}/test_result_neural_net.csv")


def test_streptomycin_called_variants_into_random_forest(capfd):
    random_forest_container = (
        "julibeg/tb-ml-random-forest-from-variants-streptomycin:v0.4.0"
    )
    variant_calling_container = "julibeg/tb-ml-variants-from-aligned-reads:v0.4.0"
    cont_arg_list = [
        [
            random_forest_container,
            "--get-target-vars -o target_vars.csv".split(),
        ],
        [
            variant_calling_container,
            f"-b {test_reads} -t target_vars.csv -o genotypes.csv".split(),
        ],
        [
            random_forest_container,
            ["genotypes.csv"],
        ],
    ]
    run_test(capfd, cont_arg_list, f"{test_data}/test_result_SM_random_forest.csv")
