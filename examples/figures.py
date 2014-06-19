#!/usr/bin/env python

__doc__ = """
This file generates the figures in the paper. In order to run this file, you
will need to first run the "analysis.py" script in this folder. This will
generate backup files in the folder sim_backup, that are used to create the
respective figures. All figures are stored as PDF and EPS files in a subfolder
"figures".


Copyright (C) 2014 Ingo Fründ

This code reproduces the analyses in the paper

    Fründ, Wichmann, Macke (2014): Quantifying the effect of inter-trial dependence on perceptual decisions. J Vis, 14(7): 9.


    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    If you use the Software for your own research, cite the paper.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import sys

if len(sys.argv)>1:
    observer = sys.argv[1]
else:
    observer = 'pk'

textfile = open ( 'results.txt', 'w' )
fullwidth=6.3
halfwidth=0.45*fullwidth
laby = 1.01

# figname = "monkeyB"
# figname = "pk"
# figname = 'gbh'
figname = 'KP'
# figname = "figure"
# figname = 'sim_nh'
# figname = 'sim'

thlev = .85

####################

import extended as tckr
tckr = tckr.Extended
from matplotlib.ticker import MaxNLocator

import threshold,history
import pylab as pl
import graphics,convenience,statistics,model
import os,cPickle,logging
import matplotlib
import limits

pl.rcParams['text.latex.preamble'] = r"""
\usepackage{helvet}
\renewcommand{\familydefault}{cmss}
\setlength{\parindent}{0em}
"""
pl.rcParams['text.usetex'] = True
pl.rcParams['legend.fontsize'] = 'small'

backup_file = os.path.join ( "sim_backup",observer+".pcl" )
results = cPickle.load ( open ( backup_file, 'r' ) )

M = results['model_nohist']
print M.applythreshold
print M.nu
print M.w
print results['model_w_hist'].w

prcformatter = matplotlib.ticker.FuncFormatter ( lambda x, pos: "%g" % (100*x,) )
myformatter = matplotlib.ticker.FuncFormatter ( lambda x, pos: "%g" % ( x, ) )

all_observers = ['CBW','ANK','CM','KKB','KP',
        'cn','fb','ip','pk','xh',
        'dc','fe','fj','gbh','kp','rz',
        'aAB','aAJ','aJR','aKP','aPK','aSS' ]
Z = 300
c_plaid,c_lum,c_jw,c_audio = [(float(x)/Z,float(y)/Z,float(z)/Z)
        for x,y,z in [(25,25,112),(34,139,34),(255,215,0),(240,255,255)]]
all_colors = [c_plaid]*5 + [c_lum]*5 + [c_jw]*6 + [c_audio]*6
all_labels = [('Plaid masking',c_plaid),('Luminance',c_lum),('Discrimination',c_jw),('Audio',c_audio)]
ecolors = list ( all_colors )
# ecolors[all_observers.index('gbh')] = 'k'
try:
    ecolors[all_observers.index(observer)] = 'y'
except:
    pass

def load_observer ( observer ):
    if observer in ["plaid","CBW","ANK","CM","KKB","KP"]:
        if observer is 'plaid':
            observer = 'CBW'

        logging.info ( "Loading data from plaid masking experiment: %s" % (observer,) )
        data,w0,plotinfo = convenience.load_plaid_c ( observer )

    elif observer in ["cn","fb","ip","pk","xh"]:

        logging.info ( "Loading data from Adelson's checker board experiment -- observer: %s" % (observer,) )
        data,w0,plotinfo = convenience.load_checker ( observer )

    elif observer in ['dc','fe','fj','gbh','kp','rz']:

        logging.info ( "Loading data from 2afc detection experiment" )
        data,w0,plotinfo = convenience.load_detection ( observer )

    elif observer[:3] == "sim":

        logging.info ( "Loading simulated data -- code: %s" %( observer[4:], ) )
        data,w0,plotinfo = convenience.load_sim ( observer[4:] )

    elif observer[0] == "a":

        observer = observer[1:]
        logging.info ( "Loading data from audio experiment, observer: %s" % \
                (observer,) )
        data,w0,plotinfo = convenience.load_audio ( observer )

    elif observer in ["nico0","nico1","nico2","konrad0","konrad1","konrad2"]:

        condition = int ( observer[-1] )
        animal    = observer[:-1]
        logging.info ( "Loading monkey %s, condition %d" % (animal, condition) )
        data,w0,plotinfo = convenience.load_monkey_condition ( animal )
        data.set_condition ( condition )

    elif observer in ["nico", "konrad"]:

        animal = observer
        logging.info ( "Loading monkey %s, all conditions" % (animal,) )
        data,w0,plotinfo = convenience.load_monkey ( animal )


    return data,w0,plotinfo

data,w0,plotinfo = load_observer ( observer )

def place_axes ( fig, x,y, fieldwidths, fieldheights, showy, margins, wh, twiny=False ):
    w,h = wh
    fieldwidths = pl.array ( fieldwidths, 'd' )
    fieldheights = pl.array ( fieldheights, 'd' )
    cwidths = pl.concatenate ( [[0],pl.cumsum ( fieldwidths )] )

    ax = []
    for i in xrange (len(fieldwidths)):
        xpos = float(x+cwidths[i]+margins[i])/w
        ypos = (y+1.)/h
        width = float(fieldwidths[i]-2*margins[i])/w
        heigh = (fieldheights[i]-2)/h
        if heigh>0:
            yax = showy[i]
            if len ( ax )>0 and twiny:
                ax.append (
                    graphics.prepare_axes (
                        fig.add_axes ( [xpos,ypos,width,heigh],
                            sharey=ax[0] ),
                        haveon=('bottom','left') if yax else ('bottom',) )
                    )
            else:
                ax.append (
                    graphics.prepare_axes (
                        fig.add_axes ( [xpos,ypos,width,heigh] ),
                        haveon=('bottom','left') if yax else ('bottom',) )
                    )
    return ax

def figure1 ( ):
    w,h = 28,15
    fig = pl.figure ( figsize=(fullwidth,h*fullwidth/w) )

    a_,b_,d_ = place_axes ( fig, 1.5,9, [9,8,9], [6]*3,
            [True]+[False]+[True], [.5,1.5,1.], (w,h) )
    c1,c2,e1,e2,e3 = place_axes ( fig, 1.5,2, [5.2,5.2,5.2,5.2,5.2], [6]*5,
            [False]*5, [.5]*5, (w,h) )
    a_.text ( .05, laby, r"\textbf{a}", transform=a_.transAxes )
    b_.text ( .05, laby, r"\textbf{b}", transform=b_.transAxes )
    d_.text ( .05, laby, r"\textbf{c}", transform=d_.transAxes )
    c1.text ( .05, laby, r"\textbf{d}", transform=c1.transAxes )
    c2.text ( .05, laby, r"\textbf{e}", transform=c2.transAxes )
    e1.text ( .05, laby, r"\textbf{f}", transform=e1.transAxes )
    e2.text ( .05, laby, r"\textbf{g}", transform=e2.transAxes )
    e3.text ( .05, laby, r"\textbf{h}", transform=e3.transAxes )


    # Figure 1A ==============================
    print "nu0",results['model_nohist'].w,results['model_nohist'].nu
    print results['model_nohist'].applythreshold
    plotinfo_ = dict ( plotinfo )
    if observer in ["KP",'sim_KP','sim_KP_nh']:
        plotinfo_['conditions'] = (0,5)
    plotinfo_['xmin'] = -25
    plotinfo_['xmax'] =  25
    if observer in ['nico','konrad']:
        plotinfo_['xmin'] = -100
        plotinfo_['xmax'] =  100
    print plotinfo_
    convenience.pmfplot ( data, results, plotinfo_, a_, errors=False )
    if observer is 'pk':
        a_.set_xlim ( -30, 30  )
        a_.set_xticks ( (-30,-20,-10,0,10,20,30) )
    a_.set_title ('')
    a_.set_ylim ( -.03, 1.03 )
    a_.yaxis.set_major_formatter ( myformatter )
    a_.xaxis.set_major_formatter ( myformatter )
    a_.set_ylabel ( "Probability" )
    a_.set_xlabel ( "Transduced stimulus intensity" )

    # Figure 1B ===============================
    textfile.write ("Figure 1B:\n" )
    l_obs,c95,caic,cpe = convenience.permutationplot ( data, results, plotinfo, b_, noaic=True )
    b_.set_title ( '' )
    b_.set_xlabel ( "Log-likelihood" )
    b_.xaxis.set_major_formatter ( myformatter )
    b_.set_xlim ( trimmed_hlim ( results['permutation_wh'][:,0], l_obs, (-2000,1000),20 ) )
    for l in b_.get_children():
        if isinstance ( l, matplotlib.legend.Legend ):
            pl.setp(l, visible=False )
    print l_obs,c95,caic,cpe,caic-6
    textfile.write ( "  l_obs = %g\n  l_95%% = %g\n  l_AIC = %g\n  cpe = %g\n" % (l_obs,c95,caic,cpe) )

    if getattr ( data, 'audio', False ):
        easy,difficult = data.performance_filter ()
        M = statistics.EvaluationCollector ( results['model_w_hist'], easy, difficult )
    else:
        M = statistics.EvaluationCollector ( results['model_w_hist'] )
    M(results['model_w_hist'])

    print results['permutation_wh'].shape

    # Figure 1Ca
    textfile.write ( "Figure 1C:\n" )
    # Variance explained on difficult trials
    hist,c95 = statistics.historytest ( results['permutation_wh'][:,3] )
    graphics.montecarlo_test ( M.vdifficult, hist, c95, ax=c1, labeling='other' )
    c1.set_xlabel ( "Var. explained [\%]\n difficult stimuli", fontsize=8, multialignment='center' )
    c1.xaxis.set_major_formatter ( prcformatter )
    c1.set_xlim ( 0, max(M.vdifficult*1.1, hist[1].max()) )
    textfile.write ( "  variance explained on difficult trials: %g, crit: %g, cpe: %g\n" %\
            (M.vdifficult,c95,pl.mean ( results['permutation_wh'][:,3]<M.vdifficult)) )

    # figure 1Cb
    # Variance explained on easy trials
    hist,c95 = statistics.historytest ( results['permutation_wh'][:,4] )
    graphics.montecarlo_test ( M.veasy, hist, c95, ax=c2, labeling='other' )
    c2.set_xlabel ( "Var. explained [\%]\n easy stimuli", fontsize=8, multialignment='center' )
    c2.xaxis.set_major_formatter ( prcformatter )
    c2.set_xlim ( 0, max(M.veasy*1.1, hist[1].max()) )
    textfile.write ( "  variance explained on easy trials: %g, crit: %g, cpe: %g\n" % \
            (M.veasy,c95,pl.mean( results['permutation_wh'][:,4]<M.veasy)) )
    textfile.write ( "  variance explained by stimulus on easy trials: %g\n" % \
            (M.vstimulus,) )

    Mh = results['model_w_hist']
    M0 = results['model_nohist']
    print "LL/trial=",(Mh.loglikelihood-results['permutation_wh'][:,0].mean())
    current_stimulus = pl.dot ( Mh.X[:,1:Mh.hf0], Mh.w[1:Mh.hf0] )
    history_features = pl.dot ( Mh.X[:,Mh.hf0:],  Mh.w[Mh.hf0:] )
    decision_signal  = current_stimulus + history_features
    textfile.write ( "  predicted slope reduction: %g\n" % ((1-0.25*pl.var(history_features)),) )
    textfile.write ( "  actual slope reductions (al_wh/al_nh): %s\n" %  (str((Mh.w[1:Mh.hf0]/M0.w[1:M0.hf0]).tolist()),) )


    # figure 1D
    S,C,V = M.stimuli,M.conditions,M.variance_explained
    conditions = pl.unique ( C )
    print conditions,plotinfo_
    for j,c in enumerate(plotinfo_['conditions']):
        print j,c,plotinfo['colors']
        i = C==plotinfo_['indices'][c][1]
        d_.plot ( S[i], V[i], '.', color=plotinfo['colors'][j] )
    pl.setp ( d_, xlim=(0,max(-plotinfo['xmin'],plotinfo['xmax'])), xlabel='Stimulus intensity',
            ylabel='Var. explained [\%]', ylim=(-.03,1.03) )
    d_.yaxis.set_major_formatter ( prcformatter )
    d_.xaxis.set_major_formatter ( myformatter )
    d_.set_xticks ( (0,10,20,30) )
    d_.set_xlim ( -.1, 30.1 )

    # Figure 1E
    textfile.write ( "Figure 1E:\n" )
    # prediction from history+stimulus
    pS_samples = results['permutation_wh'][:,6]
    pSH_samples = results['permutation_wh'][:,7]
    hist,c95 = statistics.historytest ( results['permutation_wh'][:,7] )
    if pl.var ( pSH_samples ) < 1e-7:
        mpSH = pl.mean ( pSH_samples )
        hist = (pl.array ( [len(pSH_samples)] ),pl.array([mpSH-1e-3,mpSH+1e-3]))
    graphics.montecarlo_test ( M.pSH, hist, c95, Caic=M.pstim, ax=e1, labeling='other' )
    e1.set_xlabel ( 'Prediction acc. [\%]\ndifficult stimuli', fontsize=8, multialignment='center' )
    e1.xaxis.set_major_formatter ( prcformatter )
    e1.set_xlim ( trimmed_hlim ( pSH_samples, (M.pSH,M.pstim)))
    print pSH_samples.mean(),pSH_samples.std(),"---"

    textfile.write ( "  prediction accuracy H+S, difficult: %g, crit: %g, cpe: %g\n" %\
            (M.pSH, c95, pl.mean(pSH_samples<M.pSH)) )
    textfile.write ( "  prediction accuracy S, difficult: %g, crit: %g, cpe: %g\n" %\
            (M.pstim, pl.prctile(pS_samples,95), pl.mean(pS_samples<M.pstim)) )


    hist,c95 = statistics.historytest ( results['permutation_wh'][:,5] )
    graphics.montecarlo_test ( M.phist, hist, c95, Caic=M.pstim, ax=e2, labeling='other' )
    e2.set_xlabel ( 'Prediction acc. [\%]\ndifficult stimuli', fontsize=8, multialignment='center' )
    e2.xaxis.set_major_formatter ( prcformatter )
    pH_samples = results['permutation_wh'][:,5]
    e2.set_xlim ( trimmed_hlim ( pH_samples, M.phist))
    textfile.write ( "  prection accuracy H, difficult: %g, crit: %g, cpe: %g\n" %\
            (M.phist,c95,pl.mean ( pH_samples<M.phist)) )

    hist,c95 = statistics.historytest ( results['permutation_wh'][:,8] )
    graphics.montecarlo_test ( M.peasy, hist, c95, ax=e3, labeling='other' )
    e3.set_xlabel ( 'Prediction acc. [\%]\neasy stimuli', fontsize=8, multialignment='center' )
    e3.xaxis.set_major_formatter ( prcformatter )
    peasy_samples = results['permutation_wh'][:,8]
    e3.set_xlim ( trimmed_hlim ( peasy_samples, M.peasy))
    textfile.write ( "  prection accuracy H, easy: %g, crit: %g, cpe: %g\n" %\
            (M.peasy,c95,pl.mean ( peasy_samples<M.peasy)) )

    # a_.xaxis.set_major_locator (
    #         tckr ( density=0.45, figure=fig, which=0 ) )
    b_.xaxis.set_major_locator (
            tckr ( density=0.2, figure=fig, which=0 ) )
    for ax in (c1,c2):
        ax.xaxis.set_major_locator (
                tckr ( density=0.3, figure=fig, which=0 ) )
        ax.set_xlim ( 0, None )
    d_.xaxis.set_major_locator ( tckr ( density=0.4, figure=fig,which=0 ) )
    d_.set_xlim ( -.2, None )
    for ax in (e1,e2,e3):
        ax.xaxis.set_major_locator (
                MaxNLocator ( 5 ) )
                # tckr ( density=0.3, figure=fig, which=0 ) )
    if observer in ['pk']:
        e1.set_xticks ( pl.array((60,62,64,66))/100. )
        e1.set_xlim ( .60,.67 )

    pl.savefig ( "figures/%s1.pdf" % (figname,) )
    pl.savefig ( "figures/%s1.eps" % (figname,) )

def figure2 ( ):
    w,h = 20,7
    fig = pl.figure ( figsize=(fullwidth,h*fullwidth/w) )

    a,b = place_axes ( fig, 1.1,1, [7,7,5], [6,6,0], [True,False],
            [.5,.5,1], (w,h), twiny=True )
    kl = convenience.kernelplot ( data, results, plotinfo, a, b, None )
    textfile.write ( "Figure 2\n" )
    M = results['model_w_hist']
    hr = data.gethistorykernel ( M.w[data.hf0:data.hf0+data.hlen], M.w[1] )
    hz = data.gethistorykernel ( M.w[data.hf0+data.hlen:],      M.w[1] )
    bootstrap = results['bootstrap']
    kernellen = (bootstrap.shape[1]-1)/2
    ci_stim,ci_resp,ci_corr,ci_inco = statistics.history_kernel_ci (
            bootstrap[:,kernellen:-2], bootstrap[:,:kernellen],
            hz, hr )

    textfile.write ( "  stim[-1] = %g in (%g,%g)\n" % (hz[0], ci_stim[0,0], ci_stim[1,0]) )
    textfile.write ( "  resp[-1] = %g in (%g,%g)\n" % (hr[0], ci_resp[0,0], ci_resp[1,0]) )
    textfile.write ( "  correct[-1] = %g in (%g,%g)\n" % (hz[0]+hr[0],   ci_corr[0,0], ci_corr[1,0]) )
    textfile.write ( "  incorrect[-1] = %g in (%g,%g)\n" % (-hz[0]+hr[0], ci_inco[0,0], ci_inco[1,0]) )

    a.text ( .05, laby, r"{\bf a}", transform=a.transAxes )
    b.text ( .05, laby, r"{\bf b}", transform=b.transAxes )
    pl.setp ( (a,b), title="", xlabel="Lag" )
    pl.setp ( b, yticks=(), ylabel="" )
    pl.setp ( a, ylabel="weight" )
    a.xaxis.set_major_formatter ( myformatter )
    b.xaxis.set_major_formatter ( myformatter )
    a.yaxis.set_major_formatter ( myformatter )

    handles = []
    annotations = []
    for l in kl:
        handles.append ( l )
        annotations.append ( l.get_label() )
    assert len(handles)==4
    pl.figlegend ( handles, annotations,
            'lower left', bbox_to_anchor=(15./w,2./h),
            numpoints=1,
            title='Weights of\nprevious...')

    # mx = abs ( pl.array([ci_inco, ci_stim, ci_resp, ci_corr]) ).max( )
    # mx = 0.5*10**pl.ceil(pl.log10(2*mx))
    # pl.setp ( (a,b), ylim = (-mx,mx) )
    a.yaxis.set_major_locator ( tckr ( density=0.5, figure=fig, which=1 ) )
    b.yaxis.set_visible ( False )

    if observer in ['sim_KP','sim_KP_nh']:
        pl.setp ( (a,b), ylim = (-1,1) )

    pl.savefig ( "figures/%s2.pdf" % (figname,) )
    pl.savefig ( "figures/%s2.eps" % (figname,) )

def figure3 ( ):
    w,h = 25,8.5
    fig = pl.figure ( figsize=(fullwidth,h*fullwidth/w) )

    # a,b,c,d = place_axes ( fig, 1.5,2, [9,9,5,5],[6]*4,
    #         [True]*2+[False]*2, [1.8,1.8,.5,.5], (w,h) )
    a,b,c = place_axes ( fig, 1.5,2, [9,9,5],[6]*3,
            [True]*2+[False], [1.8,1.8,.5], (w,h) )
    d = fig.add_axes ( [10,10,1,1] )
    a.text ( .05, laby, r"\textbf{a}", transform=a.transAxes )
    b.text ( .05, laby, r"\textbf{b}", transform=b.transAxes )
    c.text ( .05, laby, r"\textbf{c}", transform=c.transAxes )
    d.text ( .05, laby, r"\textbf{d}", transform=d.transAxes )
    M = results['model_w_hist']

    # Figures 3 A,B
    for condition in plotinfo['conditions']:
        condition = int ( condition )
        print "c",condition
        d_ = data.getsummary ( condition )
        # x = pl.mgrid[0:plotinfo['xmax']:100j]
        x = pl.mgrid[0:30:100j]
        # if len(data.th_features)>0:
        #     x = threshold.u_v ( x, results['model_w_hist'].nu )

        wfit  = results['model_w_hist'].w[plotinfo['indices'][condition]]
        w0fit = results['model_nohist'].w[plotinfo['indices'][condition]]
        pfit  = results['model_w_hist'].pi
        p0fit  = results['model_nohist'].pi
        x_ = threshold.u_v ( x, results['model_w_hist'].nu )
        x0 = threshold.u_v ( x, results['model_nohist'].nu )

        col = plotinfo['colors'][condition]
        pmf = 0.5*(pfit[1]+pfit[2]*model.logistic ( wfit[0]+wfit[1]*x_ )) + \
                0.5*(1-(pfit[1]+pfit[2]*model.logistic ( wfit[0]-wfit[1]*x_ )))
        p0f = 0.5*(p0fit[1]+p0fit[2]*model.logistic ( w0fit[0]+w0fit[1]*x0 )) + \
                0.5*(1-(p0fit[1]+p0fit[2]*model.logistic ( w0fit[0]-w0fit[1]*x0 )))
        print p0fit
        perror = (1-p0f-(1-pmf))/(1-p0f)

        a.plot ( x, pmf, color = col )
        a.plot ( x, p0f, color = col, linestyle='--' )
        b.plot ( x, pl.clip(perror,0,1e5), color = col )

    a.yaxis.set_major_formatter ( prcformatter )
    a.xaxis.set_major_formatter ( myformatter )
    a.set_xticks ( (0,10,20,30) )

    pl.setp ( (a,b), xlabel='Stimulus intensity' )
    a.set_ylabel ( 'Probability correct [\%]' )
    b.set_ylabel ( 'Error rate exp. [\%]' )
    b.set_xticks ( (0,10,20,30) )
    b.yaxis.set_major_locator ( tckr ( density=2, figure=fig, which=1 ) )
    b.yaxis.set_major_formatter ( prcformatter )
    b.xaxis.set_major_formatter ( myformatter )
    if observer in ['KP','sim_KP','sim_KP_nh']:
        b.set_ylim ( 0, .35 )
    if observer in ['pk']:
        pl.setp ( (a,b), xlim=(-.1,30.1) )

    # figure 3 C
    textfile.write ( "Figure 3C:\n" )
    z0 = 0
    C = statistics.EvaluationCollector ( M )
    ewh = C(results['model_w_hist'])
    enh = C(results['model_nohist'])
    hf0 = M.hf0
    # perm = results['permutation_wh']
    # # TODO: These indices have to be adapted to the revised collector
    # thresholds_wh = pl.array([C.get_thres ( perm[i,13+hf0:13+2*hf0], perm[i,12+hf0], perm[i,9:12], p=0.75 ) \
    #         for i in xrange ( 2000 )])
    # perm = results['permutation_nh']
    # thresholds_nh = pl.array([C.get_thres ( perm[i,13+hf0:13+2*hf0], perm[i,12+hf0], perm[i,9:12], p=0.75 ) \
    #         for i in xrange ( 2000 )])
    if thlev == .75:
        thind = 11
    elif thlev == .85:
        thind = 10+hf0
    else:
        raise ValueError

    for condition in xrange ( 1, M.hf0 ):
        s_wh = results['permutation_wh'][:,thind+condition]
        s_nh = results['permutation_nh'][:,thind+condition]
        # s_wh = thresholds_wh[:,condition]
        # s_nh = thresholds_nh[:,condition]
        s_ratio = s_wh/s_nh
        s_ratio_obs = ewh[thind+condition]/enh[thind+condition]
        # s_ratio_obs = results['model_w_hist'].w[condition]/results['model_nohist'].w[condition]
        z = (s_ratio_obs-pl.mean(s_ratio))/pl.std(s_ratio)
        cpe = pl.mean ( s_ratio < s_ratio_obs )
        ci = pl.prctile ( s_ratio, (2.5,97.5) )
        if z < z0 and ci[1]-ci[0] > 0:
            c0 = condition
            s_ratio_ = s_ratio
            s_ratio_obs_ = s_ratio_obs
            ci_ = ci
        textfile.write (
                "Condition %d\n  th75_ratio = %g\n  cpe = %g\n  percentiles of Null-Distribution: %g, %g\n" % \
                        (condition,s_ratio_obs,cpe,ci[0],ci[1]) )
    try:
        print "Using condition %d for figure 3C" % (c0,)
    except:
        c0 = 1
        s_ratio_ = s_ratio
        s_ratio_obs_ = s_ratio_obs
        ci_ = ci

    hist,bins = pl.histogram ( s_ratio_ )
    c.bar ( bins[:-1], hist, pl.diff ( bins ),
        edgecolor=graphics.histogram_color, facecolor=graphics.histogram_color )
    yrange = c.get_ylim ()
    # c.plot ( [1]*2, yrange, 'k:' )
    if s_ratio_obs<ci_[0]:
        c.plot ( [s_ratio_obs_]*2, (yrange[0],yrange[0]+0.85*(yrange[1]-yrange[0])), linewidth=2,
                color=graphics.observed_color )
        c.plot ( [s_ratio_obs_], [yrange[0]+0.95*(yrange[1]-yrange[0])], '*', color=graphics.observed_color )
    else:
        c.plot ( [s_ratio_obs_]*2, yrange, linewidth=2, color=graphics.observed_color )
    c.plot ( [ci_[0]]*2, yrange, color=graphics.C95_color )
    c.plot ( [ci_[1]]*2, yrange, color=graphics.C95_color )
    c.set_ylim ( *yrange )
    c.set_xlabel ( r'Threshold ratio' )
    c.xaxis.set_major_formatter ( myformatter )
    c.xaxis.set_major_formatter ( myformatter )
    # c.text ( .7, 0.7, r"$\frac{\theta_\mathrm{h}}{\theta_0}$",
    #         transform=c.transAxes )
    # c.set_xlim ( trimmed_hlim ( s_ratio_, s_ratio_obs_ ) )
    # c.xaxis.set_major_locator ( tckr ( density=0.4, figure=fig, which=0 ) )
    c.set_xlim ( .99, 1.01 )
    # c.xaxis.set_ticks ( (.95,1) )
    # c.set_xlim ( .85, 1.05 )
    c.xaxis.set_ticks ( (.99,1.,1.01) )


    # figure 3 D
    l_wh  = 0.5*results['permutation_wh'][:,[9,10]].sum(1)
    l_nh  = 0.5*results['permutation_nh'][:,[9,10]].sum(1)
    l_ratio = l_wh-l_nh
    l_ratio_obs = results['model_w_hist'].pi[[0,1]].sum()-results['model_nohist'].pi[[0,1]].sum()
    cpe = pl.mean ( l_ratio < l_ratio_obs )
    ci = pl.prctile ( l_ratio, (2.5,97.5) )
    textfile.write (
        "Figure 3D:\n  lapse_ratio = %g\n  cpe = %g\n  percentiles of Null-distribution: %g, %g\n  lapse_rate (w hist) = %g\n  lapse_rate (no hist) = %g\n" % \
                (l_ratio_obs,cpe,ci[0],ci[1],results['model_w_hist'].pi[[0,1]].sum(),results['model_nohist'].pi[[0,1]].sum()) )

    d = graphics.prepare_axes ( d, haveon=('bottom',) )
    # hist,bins = pl.histogram ( l_ratio )
    hist,bins = pl.histogram ( l_ratio, bins=good_lapse_bins ( l_ratio ) )
    # hist,bins = pl.histogram ( l_ratio, bins=pl.mgrid[-.0001:.0001:20j] )
    d.bar ( bins[:-1], hist, pl.diff(bins),
        edgecolor=graphics.histogram_color, facecolor=graphics.histogram_color, zorder=0 )
    yrange = d.get_ylim ()
    # d.plot ( [1]*2, yrange, 'k:' )
    if l_ratio_obs < ci[0] or l_ratio_obs > ci[1]:
        d.plot ( [l_ratio_obs]*2, [yrange[0], yrange[0]+0.85*(yrange[1]-yrange[0])],
                linewidth=2, color=graphics.observed_color)
        d.plot ( [l_ratio_obs], [yrange[0]+0.95*(yrange[1]-yrange[0])], '*', color=graphics.observed_color)
    else:
        print "lrobs",l_ratio_obs
        d.plot ( [l_ratio_obs]*2, yrange, color=graphics.observed_color, zorder=2)
    d.plot ([ci[0]]*2, yrange, color=graphics.C95_color, zorder=1 )
    d.plot ([ci[1]]*2, yrange, color=graphics.C95_color, zorder=1 )
    d.set_ylim ( yrange )

    d.set_xlabel ( r'Asymptote difference' )
    # d.text ( .7, 0.7, r"$\frac{\lambda_\mathrm{h}}{\lambda_0}$",
    #         transform=d.transAxes )
    # d.set_xlim ( trimmed_hlim ( l_ratio, l_ratio_obs, (0,5) ) )
    d.set_xlim ( -.003, .001 )
    d.xaxis.set_major_locator ( tckr ( density=0.4, figure=fig, which=0 ) )
    d.xaxis.set_ticks ( (-.002,0) )
    # d.set_xlim ( (.75, 1.25) )
    d.xaxis.set_major_formatter ( myformatter )

    a.set_ylim ( .49, 1.01 )

    pl.savefig ( "figures/%s3.pdf" % ( figname, ) )
    pl.savefig ( "figures/%s3.eps" % ( figname, ) )

def figure4 ( ):
    w,h = 27,9
    fig = pl.figure ( figsize=(fullwidth,h*fullwidth/w) )
    a,b,c,d,e = place_axes ( fig, 2, 2.5,
            [4]*2+[1,6,4,1,4],
            [6]*2+[0]+[6]*2+[0]+[6],
            [False]*3+[True]*4,
            [1]*7, (w,h) )
    # a,b,c = place_axes ( fig, 2, 8,   [7,7,7],[6,6,6], [False,False,True], [1.5]*3, (w,h) )
    # d,e = place_axes ( fig, 2, 1.5, [7,7,4],[6,6,0], [True]*3, [1.5]*3, (w,h) )
    a.text ( .05, laby, r"\textbf{a}", transform=a.transAxes )
    b.text ( .05, laby, r"\textbf{b}", transform=b.transAxes )
    c.text ( .05, laby, r"\textbf{c}", transform=c.transAxes )
    d.text ( .05, laby, r"\textbf{d}", transform=d.transAxes )
    e.text ( .05, laby, r"\textbf{e}", transform=e.transAxes )

    textfile.write ( "Figure 4:\n=========\n" )

    ecolors = list ( all_colors )
    # ecolors[all_observers.index('gbh')] = 'k'
    ecolors[all_observers.index(observer)] = 'y'

    if os.path.exists ( 'sim_backup/fig4_data.pcl' ):
        print "Loading fits"
        data,logstring = cPickle.load ( open ('sim_backup/fig4_data.pcl', 'r' ) )
    else:
        data = []
        logstring = []
        ldiff = []
        vhist = []
        vstim = []
        pH = []
        pSH = []
        lm = []
        th = []
        # logl_h,logl_0,vdiff,veasy,pH,pS,peasy,lm,th
        ind = [0,1,3,4,5,6,8,9,12,14]
        import scipy.io
        for o in all_observers:
            backup_file = os.path.join ( "sim_backup",o+".pcl" )
            results =  cPickle.load ( open ( backup_file, 'r' ) )
            d_ = load_observer ( o )[0]
            M = results['model_w_hist']
            if getattr ( d_, 'audio', False ):
                easy,difficult = d_.performance_filter ()
                C = statistics.EvaluationCollector ( M, easy, difficult )
            else:
                C = statistics.EvaluationCollector ( M )
            perm = results['permutation_wh']

            print "nu",M.nu,results['model_nohist'].nu

            ev = C(M)
            vstim.append ( C.vstimulus )
            vhist.append ( (C.vdifficult,pl.mean(C.vdifficult>perm[:,3])) )
            pH.append ( (C.phist,pl.mean(C.phist>perm[:,5])) )
            pSH.append ( C.pSH )
            print vhist[-1]

            ev_ = C(results['model_nohist'])
            N = len(M.r)
            logstring.append ( "  %s: l_obs=%g, l_aic=%g, cpe=%g, N=%d" % (o,ev[0],ev_[0]+6,pl.mean(ev[0]>perm[:,0]),N) )
            dat = ev[ind]
            ldiff.append ( (ev[0]-pl.mean(perm[:,0]))/N )
            dat[0] /= N
            dat[1] = (ev_[0]+6)/N
            # dat[1] = pl.prctile ( results['permutation_wh'][:,0], (95,) )/N
            # dat[1] = pl.mean(perm[:,0])/N
            lm_ = ev[[9,10]].sum()-ev_[[9,10]].sum()
            print "lm",ev[[9,10]].sum(),ev_[[9,10]].sum(),lm_
            lm.append ( (lm_,pl.mean(perm[:,9]-results['permutation_nh'][:,9]<lm_)) )
            dat[7] = lm_
            print dat[-2]

            # find optimal
            z0 = 100
            za = -100
            C = statistics.EvaluationCollector ( M )
            ewh = C(results['model_w_hist'])
            enh = C(results['model_nohist'])
            hf0 = C.hf0
            if thlev == .75:
                thind = 11
            elif thlev == .85:
                thind = 10+hf0
            else:
                raise ValueError
            aind = 13+2*(hf0-1)
            print "aind",aind
            for condition in xrange ( 1, M.hf0 ):
                s_wh = results['permutation_wh'][:,thind+condition]
                s_nh = results['permutation_nh'][:,thind+condition]
                s_ratio = s_wh/s_nh
                s_ratio_obs = ewh[thind+condition]/enh[thind+condition]
                # s_ratio_obs = results['model_w_hist'].w[condition]/results['model_nohist'].w[condition]
                z = (s_ratio_obs-pl.mean(s_ratio))/pl.std(s_ratio)
                cpe = pl.mean ( s_ratio < s_ratio_obs )
                ci = pl.prctile ( s_ratio, (2.5,97.5) )
                if s_ratio_obs < z0 and ci[1]-ci[0] > 0:
                    c0 = condition
                    s_ratio_ = s_ratio
                    s_ratio_obs_ = s_ratio_obs
                    ci_ = ci
                    z0 = s_ratio_obs

                a_wh = results['permutation_wh'][:,aind+condition]
                a_nh = results['permutation_nh'][:,aind+condition]
                a_ratio = a_wh/a_nh
                a_ratio_obs = ewh[aind+condition]/enh[aind+condition]
                print ewh[aind+condition],enh[aind+condition],
                z = (a_ratio_obs-pl.mean(a_ratio))/pl.std(a_ratio)
                print z
                cpe = pl.mean ( a_ratio < a_ratio_obs )
                ci = pl.prctile ( a_ratio, (2.5,97.5) )
                if a_ratio_obs>za and ci[1]-ci[0] > 0:
                    a0 = condition
                    a_ratio_ = a_ratio
                    a_ratio_obs_ = a_ratio_obs
                    cia = ci
                    za = a_ratio_obs

            dat[8] = s_ratio_obs_
            dat[9] = a_ratio_obs_
            th.append ( (s_ratio_obs_, pl.mean ( s_ratio_<s_ratio_obs ) ) )

            data.append ( dat )
        data = pl.array ( data )

        mdict = {'prediction_acc_from_full_model': pSH}
        scipy.io.savemat ( 'prediction_full_model.mat', mdict )

        vhist = pl.array ( vhist )
        pH = pl.array ( pH )
        lm = pl.array ( lm )
        print "lm = ", repr(lm)
        print "a =", repr(data[:,9])
        th = pl.array ( th )

        logstring.append ( "  average increase in log_likelihood per trial: %g +/- %g" %(pl.mean (ldiff), pl.std(ldiff)) )
        logstring.append ( "  average increase in log_likelihood beyond AIC: %g +/- %g" % \
                (pl.mean(data[:,0]-data[:,1]),pl.std(data[:,0]-data[:,1])) )
        logstring.append ( "  variance explained by history: %g +/- %g , (significant: %d/%d, max: %g)" %\
                (pl.mean(vhist[:,0]),pl.std(vhist[:,0]),pl.sum(vhist[:,1]<.05),vhist.shape[0], vhist[:,0].max() ) )
        logstring.append ( "  variance explained by stimulus: %g +/- %g" %(pl.mean(vstim),pl.std(vstim)) )
        logstring.append ( "  prediction pS: %g +/- %g" % (pl.mean(data[:,5]),pl.std(data[:,5])) )
        logstring.append ( "  prediction pH: %g +/- %g (%d/%d significant)" % \
                (pl.mean(data[:,4]),pl.std(data[:,4]),pl.sum(pH[:,1]<0.05),pH.shape[0]) )
        logstring.append ( "  prediction pSH: %g +/- %g" % (pl.mean(pSH),pl.std(pSH)) )
        logstring.append ( "  prediction on easy trials: %g +/- %g" % (pl.mean(data[:,6]),pl.std(data[:,6])) )
        logstring.append ( "  variance on easy trials: %g +/- %g" % (pl.mean(data[:,3]),pl.std(data[:,3])) )

        logstring.append ( "  asymptote ratio: %g +/- %g (significant: %d/%d)" %\
                (pl.mean(lm[:,0]),pl.std(lm[:,0]),pl.sum(lm[:,1]<0.05), lm.shape[0] ) )
        logstring.append ( "  threshold ratio: %g +/- %g (significant: %d/%d)" %\
                (pl.mean(th[:,0]),pl.std(th[:,0]),pl.sum(th[:,1]<0.05), th.shape[0] ) )
        print th

        logstring = "\n".join ( logstring )

        cPickle.dump ( (data,logstring), open ( 'sim_backup/fig4_data.pcl', 'w' ) )
    textfile.write ( logstring )
    data = pl.array ( data )

    # a.scatter ( data[:,0], data[:,1], c=all_colors, edgecolor=ecolors )
    boxdotplot ( a, data[:,0]-data[:,1], color=all_colors, ecolor=ecolors, jitter=0.05 )
    for o,d_ in zip ( all_observers, data[:,0]-data[:,1] ):
        print o,d_
    a.set_xlim ( -.5, 2 )
    # a.plot ( [0,-1],[0,-1], 'k:' )
    a.set_ylabel ( "Likelihood increase\n from history", multialignment='center', fontsize=8 )
    # a.set_ylabel ( r'95th percentile likelihood', horizontalalignment='center', fontsize=8 )
    # a.set_xlim ( a.set_ylim ( -.7, -.2 ) )
    a.yaxis.set_major_locator ( tckr ( density=0.4, figure=fig, which=0 ) )
    a.yaxis.set_major_formatter ( myformatter )
    # a.set_xticks ( (0,.1,.2) )

    boxdotplot ( b, data[:,2], color=all_colors, ecolor=ecolors, jitter=0.05 )
    b.set_ylabel ( 'Var. explained [\%]', fontsize=8 )
    b.yaxis.set_major_formatter ( prcformatter )
    b.set_xlim ( -.5, 2)
    b.set_ylim ( -.03, .53 )
    textfile.write ( "\nFigure 4b:\n" )
    textfile.write ( "  variance explained by history on difficult trials: %g +/- %g\n" % ( pl.mean ( data[:,2] ), pl.std(data[:,2]) ) )

    c.scatter ( data[:,4], data[:,5], s=5, c=all_colors, edgecolor=ecolors )
    c.plot ( [0,1],[0,1], 'k:' )
    c.set_xlabel ( 'Prediction acc.\n from history',fontsize=8, multialignment='center' )
    c.set_ylabel ( 'Prediction acc.\n from stimulus', multialignment='center',fontsize=8 )
    c.set_xlim ( c.set_ylim ( 0.5, .7 ) )
    c.xaxis.set_major_formatter ( prcformatter )
    c.yaxis.set_major_formatter ( prcformatter )


    # d.plot ( [0,1.2],[1,1],  'k:', zorder=0 )
    # d.plot ( [1,1], [0,1.2], 'k:', zorder=0 )
    # d.scatter ( data[:,7], data[:,8], c=all_colors, edgecolor=ecolors )
    boxdotplot ( d, data[:,7], color=all_colors, ecolor=ecolors, jitter=0.05 )
    boxdotplot ( e, data[:,8], color=all_colors, ecolor=ecolors, jitter=0.05 )
    # d.set_xlabel ( 'Asymptote ratio',fontsize=8 )
    # d.set_ylabel ( 'Threshold ratio', horizontalalignment='center',fontsize=8 )
    d.set_ylabel ( "Asymptote difference" )
    e.set_ylabel ( "Threshold ratio" )

    print "lm =",repr(data[:,7])

    print "THR=",repr(data[:,8])
    textfile.write ( "  max threshold ratio: %g\n" % (data[:,-1].max(),) )
    # d.xaxis.set_major_locator ( tckr ( density=.5, figure=fig, which=0 ) )
    # d.yaxis.set_major_locator ( tckr ( density=.5, figure=fig, which=1 ) )
    d.set_ylim ( -.05, .05 )
    # e.set_ylim ( 0, 2 )
    # d.xaxis.set_major_formatter ( myformatter )
    d.yaxis.set_major_formatter ( myformatter)
    e.yaxis.set_major_formatter ( myformatter )

    p = []
    l = []
    dummy_axes = fig.add_axes ( [2,2,1,1], xticks=(), yticks=() )
    for e,c in all_labels:
        p.append ( dummy_axes.plot( [0],[0],'o', color=c, markeredgecolor=c ) )
        l.append ( e )
    # igbh = all_observers.index ( 'gbh' )
    iobs = all_observers.index ( observer )
    # p.append ( dummy_axes.plot ( [0],[0], 'o', color=all_colors[igbh], markeredgecolor=ecolors[igbh] ) )
    # l.append ( 'J\&W (2006), gbh' )
    p.append ( dummy_axes.plot ( [0],[0], 'o', color=all_colors[iobs], markeredgecolor=ecolors[iobs] ) )
    l.append ( 'Observer KP' )

    fig.legend ( p, l, bbox_transform=fig.transFigure, numpoints=1, ncol=5, loc=8 )

    fig.savefig ( 'figures/%s4.pdf' % (figname,) )
    fig.savefig ( 'figures/%s4.eps' % (figname,) )

def figure5 ():
    w,h = 29,8
    fig = pl.figure ( figsize=(fullwidth,h*fullwidth/w) )
    a,b,c = place_axes ( fig, 2, 1, [8,1,6,1,6,6], [6,0,6,0,6,0], [True]*6, [1.]*6, (w,h) )
    dummy_axes = fig.add_axes ( [2,2,1,1], xticks=(), yticks=() )
    a.text ( .05, laby, r"\textbf{a}", transform=a.transAxes )
    b.text ( .05, laby, r"\textbf{b}", transform=b.transAxes )
    c.text ( .05, laby, r"\textbf{c}", transform=c.transAxes )

    if os.path.exists ( 'sim_backup/all_kernels.pcl' ):
        print "Loading kernels"
        kernels = cPickle.load ( open ('sim_backup/all_kernels.pcl', 'r' ) )
    else:
        hf = data.h
        kernels = []
        for o in all_observers:
            backup_file = os.path.join ( "sim_backup",o+".pcl" )
            results_ =  cPickle.load ( open ( backup_file, 'r' ) )
            M = results_['model_w_hist']
            C = statistics.Kernel_and_Slope_Collector (
                    hf, M.hf0, range(1,M.hf0) )
            kernels.append ( C(M) )
        kernels = pl.array ( kernels ).T
        cPickle.dump ( kernels, open ( 'sim_backup/all_kernels.pcl', 'w' ) )

    kl = convenience.kernelplot ( data, results, plotinfo, a, dummy_axes, 'upper right' )
    a.set_title ( '' )
    a.get_legend().get_frame().set_linewidth(0)
    pl.setp ( kl, markersize=2 )
    for o in ["cn","fb","ip","pk","xh"]:
        if not o is observer:
            i = all_observers.index(o)
            al = kernels[14,i]
            kr = kernels[:7,i]*al
            kz = kernels[7:14,i]*al
            a.plot ( pl.arange ( 1, 8 ), kz, '-', color=graphics.stimulus_color )
            a.plot ( pl.arange ( 1, 8 ), kr, '-', color=graphics.response_color )
    a.set_xlabel ( 'Lag' )
    a.set_ylabel ( "Equivalent contrast [\%]" )
    a.set_ylim ( -1, 1 )
    a.set_yticks ( (-1,0,1) )
    a.set_xlim ( .5,7.5 )
    a.xaxis.set_major_formatter ( myformatter )
    a.yaxis.set_major_formatter ( myformatter )

    b.scatter ( kernels[0,:], kernels[7,:], c=all_colors, edgecolor=ecolors )

    for i,o in enumerate ( all_observers ):
        print o,kernels[0,i], kernels[7,i]

    b.set_xlabel ( "Previous response" )
    b.set_ylabel ( "Previous stimulus" )
    b.set_xlim ( -2,1 )
    b.set_ylim ( -1,2 )
    b.set_xticks ( (-2,-1,0,1) )
    b.set_yticks ( (-1,0,1,2) )
    b.xaxis.set_major_formatter ( myformatter )
    b.yaxis.set_major_formatter ( myformatter )

    # b.set_yticklabels ( (-2,"",0,"",2) )
    # b.set_xticklabels ( (-2,"",0,"",2) )
    # b.xaxis.set_major_locator ( tckr ( density=.3, figure=fig, which=0 ) )
    # b.yaxis.set_major_locator ( tckr ( density=.3, figure=fig, which=1 ) )

    c.scatter ( kernels[1:7,:].sum(0), kernels[8:14,:].sum(0), c=all_colors, edgecolor=ecolors )
    c.set_xlabel ( "Summed response weight" )
    c.set_ylabel ( "Summed stimulus weight" )
    c.set_xlim ( c.set_ylim ( -.25,.25 ) )
    c.set_xticks ( (-1,-.5,0,.5,1) )
    c.set_yticks ( (-1,-.5,0,.5,1) )
    c.set_yticklabels ( (-1,"",0,"",1) )
    c.set_xticklabels ( (-1,"",0,"",1) )

    p = []
    l = []
    for e,c in all_labels:
        p.append ( dummy_axes.plot( [0],[0],'o', color=c, markeredgecolor=c ) )
        l.append ( e )
    iobs = all_observers.index ( observer )
    p.append ( dummy_axes.plot ( [0],[0], 'o', color=all_colors[iobs], markeredgecolor=ecolors[iobs] ) )
    l.append ( 'Observer KP' )
    leg = fig.legend ( p, l, loc='right', bbox_transform=fig.transFigure, numpoints=1 )
    leg.get_frame().set_linewidth(0)

    fig.savefig ( 'figures/%s5.pdf' % (figname,) )
    fig.savefig ( 'figures/%s5.eps' % (figname,) )

def figure5_ ():
    fig = pl.figure ( figsize=(halfwidth,halfwidth) )

    if os.path.exists ( 'sim_backup/all_kernels.pcl' ):
        print "Loading kernels"
        kernels = cPickle.load ( open ('sim_backup/all_kernels.pcl', 'r' ) )
    else:
        kernels = []
        for o in all_observers:
            backup_file = os.path.join ( "sim_backup",o+".pcl" )
            results =  cPickle.load ( open ( backup_file, 'r' ) )
            M = results['model_w_hist']
            se = results['bootstrap'][:,:14].std(0).mean()
            C = statistics.KernelCollector ( history.history_impulses(), M.hf0, 1 )
            kernels.append ( C(M)/se )
        kernels = pl.array ( kernels ).T
        cPickle.dump ( kernels, open ( 'sim_backup/all_kernels.pcl', 'w' ) )

    a = graphics.prepare_axes ( fig.add_axes ( [.2, .2, .7, .7] ) )

    print kernels.shape

    a.scatter ( kernels[0,:], kernels[7,:], c=all_colors, edgecolor=all_colors )

    a.fill ( [-100,0,0,-100],[0,0,100,100],
            facecolor=[1.,.9,.9],
            edgecolor=[1.,.9,.9], zorder=0 )
    a.fill ( [100,0,0,100],[0,0,-100,-100],
            facecolor=[.9,.9,1.],
            edgecolor=[.9,.9,1.], zorder=0 )

    a.set_xlabel ( "previous response" )
    a.set_ylabel ( "previous stimulus" )

def figure6 ():
    w,h = 22,7
    monkeys = ["nico2","konrad2"]
    monkeyresults = [ cPickle.load ( open ( 'sim_back/'+monkey+'.pcl' ) ) for monkey in monkeys ]


def plot_allkernels ( hr, hz, ax ):
    ax.plot ( hr, 'o:', color=graphics.response_color,
            markeredgecolor=0.7*graphics.response_color )
    ax.plot ( hz, 'o:', color=graphics.stimulus_color,
            markeredgecolor=0.7*graphics.stimulus_color )
    ax.plot ( -hz+hr, 'o:', color=graphics.incorrect_color,
            markeredgecolor=0.7*graphics.incorrect_color )
    ax.plot ( hz+hr, 'o:', color=graphics.correct_color,
            markeredgecolor=0.7*graphics.correct_color )


def set_ax_hlims ( ax, dist, obs ):
    e = dist.ptp()
    dlims = [dist.min(),dist.max()]
    ticks = pl.array(dlims + [obs])
    lims = min(ticks)-.05*e,max(ticks)+.05*e
    ax.set_xlim ( *lims )
    ax.set_xticks ( ticks.round(2) )

def trimmed_hlim ( dist, point, rng=(0,1), step=0.05 ):
    hmin,hmax = rng
    while pl.mean ( dist<hmin ) <.01 and (pl.array(point)>hmin).all():
        hmin += step
    while pl.mean ( dist>hmax ) <.01 and (pl.array(point)<hmax).all():
        hmax -= step
    return min(hmin - step,pl.array(point).min()-step), max(hmax + step,pl.array(point).max()+step)

def good_lapse_bins ( dist ):
    hmin,hmax = trimmed_hlim ( dist, point=dist.mean(), rng=(0,5), step=0.001 )
    return pl.mgrid[hmin:hmax:20j]

def boxdotplot ( ax, data, color, ecolor, jitter=0.01, yshift=0 ):
    yc = jitter*pl.randn ( len(data) ) + yshift
    ax.scatter ( yc, data, s=5, c=color, edgecolor=ecolor )
    prc = pl.prctile ( data )
    ax.plot ( [yshift+1]*2, [prc[0],prc[-1]], 'k-', zorder=0 )
    jitter *= 5
    ax.fill ( yshift+pl.array([1-jitter,1+jitter,1+jitter,1-jitter]),
            [prc[1],prc[1],prc[3],prc[3]],
            facecolor='w', edgecolor='k', zorder=1 )
    ax.plot ( yshift+pl.array([1-jitter,1+jitter]), [prc[2]]*2, 'k-', zorder=2 )
    ax.plot ( [yshift+1], [pl.mean(data)], 'o', color='w', markeredgecolor='k' )

    ax = graphics.prepare_axes ( ax, haveon=('left',) )

if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    figure4 ()
    figure5 ()
    figure6 ()
    textfile.close()
    pl.show()
