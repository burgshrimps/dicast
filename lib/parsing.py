import argparse
import sys

def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='dicast')

    subparsers = parser.add_subparsers(dest='command', help='Subcommands')

    parser_prepare = subparsers.add_parser('prepare', help='Variant preparation and feature extraction')
    parser_prepare.add_argument('sample', help='Sample name')
    parser_prepare.add_argument('ref', help='Reference genome name', default='hg38')
    parser_prepare.add_argument('params', help='Preparation parameter file')
    parser_prepare.add_argument('workdir', help='Working and output directory')

    parser_train = subparsers.add_parser('train', help='Train a new model')
    parser_train.add_argument('svtype', help='SV type to train the model for')
    parser_train.add_argument('params', help='Training parameter file')
    parser_train.add_argument('clfparams', help='Model parameter file')
    parser_train.add_argument('clfname', help='Name of model to train')
    parser_train.add_argument('--chr_excl', help='Chromsomes to exclude from training', nargs='+', default=[])

    parser_test = subparsers.add_parser('test', help='Test a model')
    parser_test.add_argument('svtype', help='SV type to test the model for')
    parser_test.add_argument('clf', help='Input file to load the model from')
    parser_test.add_argument('params', help='Test parameter file')
    parser_test.add_argument('workdir', help='Working and output directory')
    parser_test.add_argument('--chr_incl', help='Chromsomes to include for testing', nargs='+', default=['all'])

    parser_predict = subparsers.add_parser('predict', help='Make predictions for a set of variants')
    parser_predict.add_argument('svtype', help='SV type to make predictions for')
    parser_predict.add_argument('clf', help='Input file to load the model from')
    parser_predict.add_argument('params', help='Prediction parameter file')
    parser_predict.add_argument('output', help='Output file to save the predictions')

    parser_curate = subparsers.add_parser('curate', help='Generate a set of variants for manual curation')
    parser_curate.add_argument('svtype', help='SV type to curate')
    parser_curate.add_argument('params', help='Curation parameter file')
    parser_curate.add_argument('clfparams', help='Model parameter file')
    parser_curate.add_argument('clfname', help='Name of model to train')

    return parser.parse_args(arguments)