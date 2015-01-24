#!/usr/local/epd/bin/python

__author__ = "Brandon Ayers"
__copyright__ = "Copyright 2013, Shantanu H. Joshi, Brandon Ayers, \
                 Ahmanson-Lovelace Brain Mapping Center, University of California Los Angeles"
__email__ = "ayersb@ucla.edu"

import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Compute thickness and register batch of top and bottom callosal curves.\n")
    parser.add_argument('subjectIDs', help='Ordered list of subject identifiers')
    parser.add_argument('topCurves', help='Ordered list of paths to top callosal segmentation files')
    parser.add_argument('botCurves', help='Ordered list of paths to bottom callosal segmenation files')
    parser.add_argument('-odir', dest='odir', help='output directory', required=True)
    parser.add_argument('-templateID', dest='templateID', help='SubjectID associated with curves to be used as template for group registration', default=False, required=False)
    parser.add_argument('-linearTemplate', dest='linearTemplate', help='Do not use elastic matching if template id is specified.', action='store_true', default=False, required=False)
    parser.add_argument('-listInput', dest='listInput', help='Input curves names directly into arguments', action='store_true', default=False, required=False)
    parser.add_argument('-open', dest='open_curves', help='match open curves', action='store_true', default=False, required=False)
    parser.add_argument('-linear', dest='linear', help='use linear matching', action='store_true', default=False, required=False)
    parser.add_argument('-norotate', dest='norotate', help='do not align rotations', action='store_true', default=False,required=False)
    parser.add_argument('-noplot', dest='noplot', help='do not save graphical plots', action='store_true', default=False,required=False)
    parser.add_argument('-resize', dest='resize', help='resize to the specified number of vertices', required=False, default=100)
    parser.add_argument('-altThicknessReg', dest='altReg', action='store_true', help='Use alternative method to register thickness values', required=False, default=False)
    args = parser.parse_args()
    corpus_callosum_analyze(args.subjectIDs, args.topCurves, args.botCurves, args.odir, args.templateID, args.listInput,
                            args.open_curves, args.linear, args.norotate, args.noplot, args.resize,
                            args.linearTemplate, args.altReg)


def corpus_callosum_analyze(subject_ids, top_curves, bot_curves, odir, template_id, list_input=False,
                            open_curves=False, linear=False, no_rotate=False, no_plot=False, resize=300, linear_template_matching=False, altReg=False):
    """ Takes as input
    """
    from ccshape.corpus_callosum import CorpusCallosumThickness
    if not list_input:
        subject_ids = open(os.path.abspath(subject_ids))
        subject_ids = subject_ids.read().split()
        subject_ids.sort()
        top_curves = open(os.path.abspath(top_curves))
        top_curves = top_curves.read().split()
        top_curves.sort()
        bot_curves = open(os.path.abspath(bot_curves))
        bot_curves = bot_curves.read().split()
        bot_curves.sort()

    if template_id:
        if template_id not in subject_ids:
            raise ValueError("Invalid template ID!")
        if not os.path.exists(os.path.join(odir, 'template')):
            os.makedirs(os.path.join(odir, 'template'))
        template_index = subject_ids.index(template_id)
        template_cc = CorpusCallosumThickness(template_id, top_curves[template_index], bot_curves[template_index],
                                     linear=linear, outdir=os.path.join(odir, 'template'),
                                     linear_template_matching=linear_template_matching, alt_registration=altReg)

        template_cc.output_ucf()
        if not no_plot:
            template_cc.plot()
            if not linear and not linear_template_matching:
                template_cc.plot_comparison()
        template_cc.save_matching_data()

        if linear:
            template_curve = template_cc.joined_curve_linear
        else:
            template_curve = template_cc.joined_curve_elastic
    else:
        template_curve = False

    for i in xrange(len(subject_ids)):
        if template_id:
            if i == template_index:
                continue
        if not (os.path.isfile(top_curves[i]) or os.path.isfile(bot_curves[i])):
            raise ValueError(subject_ids[i]+" is missing a curve file!")
        if not os.path.exists(os.path.join(odir, subject_ids[i])):
            os.makedirs(os.path.join(odir, subject_ids[i]))
        current_cc = CorpusCallosumThickness(subject_ids[i], top_curves[i], bot_curves[i], linear=linear,
                                    template_curve=template_curve, outdir=os.path.join(odir, subject_ids[i]),
                                    linear_template_matching=linear_template_matching, alt_registration=altReg)
        current_cc.output_ucf()
        if not no_plot:
            current_cc.plot(plot_linear=linear)
            if not linear and not linear_template_matching:
                current_cc.plot_comparison()
        current_cc.save_matching_data()

if __name__ == '__main__':
    main()








