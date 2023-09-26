import argparse
import sys

from drug.__init__ import __VERSION__, ASSAY_DICT
import drug.toolkits.utils as utils


class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def main():
    """celescope cli
    """
    parser = argparse.ArgumentParser(prog='drug',
                                     description='A transcriptome sequencing analysis tool based on barcode and UMI.',
                                     formatter_class=ArgFormatter)
    parser.add_argument('-v', '--version', action='version', version=__VERSION__)
    subparsers = parser.add_subparsers(help="sub-command help")

    for assay in ASSAY_DICT:
        text = ASSAY_DICT[assay]
        subparser_1st = subparsers.add_parser(name=assay,
                                              help=text)

        # add 2ed subparser
        subparser_2ed = subparser_1st.add_subparsers(prog=assay,
                                                     description=ASSAY_DICT[assay])

        # import __STEPS__
        init_module = utils.find_assay_init(assay)
        __STEPS__ = init_module.__STEPS__

        for step in __STEPS__:
            # import function and opts
            step_module = utils.find_step_module(assay, step)
            func = getattr(step_module, step)
            func_opts = getattr(step_module, f"get_{step}_para")
            parser_step = subparser_2ed.add_parser(step, formatter_class=ArgFormatter)
            func_opts(parser_step, optional=True)
            parser_step.set_defaults(func=func)


        
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
    else:
        args.func(args)


if __name__ == '__main__':
    main()