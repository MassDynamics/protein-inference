import argparse
from protein_inference.protein_inference_runner import ProteinInferenceRunner

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Protein Inference',
        description='Infers proteins from psms'
    )

    parser.add_argument("--output-directory",
                        dest='output_directory',
                        help='Directory of where output is place',
                        required=True)
    parser.add_argument("--target-path",
                        dest='target_path',
                        help='path to target file',
                        required=True)
    parser.add_argument("--decoy-path",
                        dest='decoy_path',
                        help='path to decoy file',
                        required=True)
    parser.add_argument("--psm-q-value-threshold",
                    dest='psms_q_value_threshold',
                    help='filtering_applied_to_input_psms',
                    default=0.01,
                    required=False)

    args = parser.parse_args()

    ProteinInferenceRunner().run(args.target_path, args.decoy_path,
                                 args.output_directory,
                                 psms_q_value_threshold = float(args.psms_q_value_threshold))
