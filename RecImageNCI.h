//
//
//

#import <Foundation/Foundation.h>
#import <Accelerate/Accelerate.h>

@class RecImage;

// p-image mode
enum {
	NCI_REAL = 0,
	NCI_CPX_MEAN,
	NCI_CPX_VAR
};

// c-func
void		nciCalcColor(int mode, float x, float y, float *r, float *g, float *b);
void		Rec_est_poly_2d(float *coef, int ordx, int ordy, float *re, float *im, int xDim, int yDim);
void		dump_hist(RecImage *img, int n, char *path);
void		dump_tdist(int df);

@interface RecImage (NCI)

- (RecLoop *)blockAvgForLoop:(RecLoop *)lp blockSize:(int)blkSz dummy:(int)ds;
- (RecLoop *)breakLoop:(RecLoop *)avgLp blockSize:(int)sz dummy:(int)ds;
- (RecImage *)kFilt:(int)mmt;	// mmt 0:dipole, 1:quadrupole
- (RecImage *)toDipole;
- (RecImage *)toQuadrupole;

- (void)removeDC2dWithMask:(RecImage *)mask;

//- (RecImage *)reorderWithDda:(int)dda phsLp:(RecLoop *)phsLp avgLp:(RecLoop *)avgLp phsArray:(int *)array;
- (RecImage *)reorderWithDda:(int)dda phsLp:(RecLoop *)phsLp avgLp:(RecLoop *)avgLp phsArray:(int *)array;
- (RecImage *)removeDdaAndMakePo2:(int)dda;
- (RecImage *)addPhs:(RecLoop *)phsLp protocol:(int)protocol dda:(int)dda;	// reverse of above
- (RecImage *)phsSliceAtIndex:(int)ix nPhs:(int)nphs dda:(int)dda;
- (RecImage *)phsSliceAtIndex:(int)ix withPhsArray:(int *)phsArray dda:(int)dda;
- (void)driftCorr;
- (RecImage *)reorderBoldWithDda:(int)dda nEvent:(int)ne evLen:(int)el slcLp:(RecLoop *)slcLp phsLp:(RecLoop *)phsLp avgLp:(RecLoop *)avgLp delay:(int)del;
- (RecImage *)makeSubSeriesWithDda:(int)dda phsLp:(RecLoop *)phsLp stim:(int)ix1 ref:(int)ix2;
- (void) addTestTag:(int)dda nPhs:(int)nPhs phsArray:(int *)phsArray amp:(float)amp;	// put test marker
- (RecImage *)filtImag:(float)w forLoop:(RecLoop *)lp;	// HPfilt imag (frac = 1), avg real
- (void)clearSliceAtIndex:(int)ix;
- (void)pcaFilt:(int)ns;	// PCA type variance filter
- (void)pcaFilt:(int)ns stim:(float *)st;	// PCA type variance filter
- (void)outerProd:(RecImage *)wt andVector:(float *)bs len:(int)len;	// ## make basis RecImage
- (void)singleFreqFilt:(int)freq forLoop:(RecLoop *)lp;	// remove single frequency (single fourier coeff)
- (void)cycFilt:(int)len forLoop:(RecLoop *)lp;	// remove single frequency (cyclic filter)
- (void)accumFilt:(int)len forLoop:(RecLoop *)lp;	// accum filter
- (void)leakyIntWithTau:(float)tau cutOff:(float)freq forLoop:(RecLoop *)lp;

- (void)sinFilt:(int)mode ang:(float)th;    // arg chenged
//- (RecImage *)cosRot:(int)nTh;
- (void)imagToPhaseUsingMag:(RecImage *)mg;
- (void)windowForLoop:(RecLoop *)lp mode:(int)filt width:(float)sd;

- (void)baseline;   // remove low freq phase
- (RecImage *)pImageWithRef:(RecImage *)ref forLoop:(RecLoop *)avgLp thres:(float)thres mode:(int)mode;
- (RecImage *)fusionWith:(RecImage *)sg gain:(float)g mode:(int)mode;
- (RecImage *)fusionWithPimg:(RecImage *)sg gain:(float)g;
- (void)timePcorr:(RecLoop *)lp;

- (void)pestPoly2d:(RecImage *)coef;
- (void)baseline2Avg:(RecLoop *)avgLp phs:(RecLoop *)phs;
//- (void)baseline3Rep:(RecLoop *)repLp;
- (void)quadFilt:(float)wt;
- (void)spikeFilt:(float)thres;
- (void)removeDC;	// 1d version
- (RecImage *)dsbDecForLoop:(RecLoop *)lp width:(float)w center:(float)c;
- (void)cos1dForLoop:(RecLoop *)lp cyc:(int)cyc power:(int)pw lowPass:(BOOL)lpf;
- (void)fCos1dForLoop:(RecLoop *)lp cyc:(int)cyc power:(int)pw lowPass:(BOOL)lpf;
- (void)gammaFiltForLoop:(RecLoop *)lp width:(float)w;

//	####
- (RecImage *)meanDiffWithRef:(RecImage *)ref forLoop:(RecLoop *)lp;
- (RecImage *)sdDiffWithRef:(RecImage *)ref forLoop:(RecLoop *)lp;

// "difference", or make diff after processing
- (RecImage *)difForLoop:(RecLoop *)phs stim:(int)ix1 ref:(int)ix2;
- (RecImage *)difForLoop:(RecLoop *)phs stim:(int)ix1 ref1:(int)ix2 ref2:(int)ix3;	// use 2 phs for higher SNR
- (void)varFiltForLoop:(RecLoop *)lp thres:(float)th;

// "differentiate" along loop
- (RecImage *)dxForLoop:(RecLoop *)lp;

// statistics
- (RecImage *)pImageWithRef:(RecImage *)ref complex:(BOOL)cpx;	// t or F, zLoop
- (RecImage *)pImageWithStim:(int)ix1 ref:(int)ix0 complex:(BOOL)cpx;	// t or F-stat
- (void)f2pDf1:(int)df1 df2:(int)df2;
- (void)t2pDf:(int)df;
- (void)pvalThres:(float)th;

// complex version of (RecImage *)varForLoop: ### varForLoop should work for cpx too ### chk
//- (RecImage *)cvarForLoop:(RecLoop *)lp;
//- (RecImage *)cvarForLoop:(RecLoop *)lp withMean:(RecImage *)mn;
- (RecImage *)phsAvgForLoop:(RecLoop *)lp;
- (RecImage *)phsVarForLoop:(RecLoop *)lp withMean:(RecImage *)mn;
- (RecImage *)phsSDForLoop:(RecLoop *)lp;
- (RecImage *)sdMaskForLoop:(RecLoop *)lp andLoop:(RecLoop *)lp2 thres:(float)thres;
- (RecImage *)sdMaskForLoop:(RecLoop *)lp thres:(float)thres;

// color ref
- (RecImage *)colMapImage:(int)mode;	// 1:diple, 2:quadrupole


// debug
- (void)dumpPhaseSDForLoop:(RecLoop *)lp image:(NSString *)path histo:(NSString *)histo;

@end

// BOLD response curve
RecImage *genResp(RecLoop *lp, float dt);
