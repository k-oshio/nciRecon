//
//	RecImage(NCI)
//

#import <RecKit/RecKit.h>
#import "RecImageNCI.h"
#import <NumKit/NumKit.h>

#import "timer_macros.h"

@implementation RecImage (NCI)

- (RecLoop *)breakLoop:(RecLoop *)avgLp blockSize:(int)nBlk dummy:(int)ds
{
    void			(^proc)(float *p, float *q, int skip, int len);
    RecImage        *img;
    RecLoop         *blkLp, *newAvgLp;
    RecLoopControl  *blkLc, *lc = [self control];
	int				blkLen;
	int				nAvg;

// create new loop structure and img
    nAvg = [avgLp dataLength];
    blkLen = (nAvg - ds) / nBlk;
    blkLp = [RecLoop loopWithName:@"blk" dataLength:nBlk];
	newAvgLp = [RecLoop loopWithName:@"newAvg" dataLength:blkLen];
	blkLc = [RecLoopControl controlWithControl:lc];
	[blkLc replaceLoop:avgLp withLoop:newAvgLp];
	[blkLc insertLoop:blkLp beforeLoop:newAvgLp];

	img = [RecImage imageOfType:RECIMAGE_COMPLEX withControl:blkLc];

	proc = ^void(float *p, float *q, int skip, int len) {
		float	*bp, *bq;
        int		j, k, srcIx, dstIx;

		bp = [img currentDataWithControl:blkLc];
		bq = bp + [img dataLength];
		for (k = 0; k < nBlk; k++) {
		// copy to new loop
			for (j = 0; j < blkLen; j++) {
				srcIx = (k * blkLen + j + ds) * skip;
				dstIx = (k * blkLen + j) * skip;
				bp[dstIx] = p[srcIx];
				bq[dstIx] = q[srcIx];
			}
		}
    };
    [self applyComplex1dProc:proc forLoop:avgLp control:lc];

    [self copyIvarOf:img];

    return newAvgLp;
}

- (RecLoop *)blockAvgForLoop:(RecLoop *)avgLp blockSize:(int)nBlk dummy:(int)ds
{
    void			(^proc)(float *p, float *q, int skip, int len);
    RecImage        *img;
    RecLoop         *blkLp;
    RecLoopControl  *blkLc, *lc = [self control];
	int				blkLen;
	int				nAvg;
	int				blkSkip;

// create new loop structure and img
    nAvg = [avgLp dataLength];
    blkLen = (nAvg - ds) / nBlk;
    blkLp = [RecLoop loopWithName:@"blk" dataLength:nBlk];
	blkLc = [RecLoopControl controlWithControl:lc];
	[blkLc replaceLoop:avgLp withLoop:blkLp];
	img = [RecImage imageOfType:RECIMAGE_COMPLEX withControl:blkLc];
	blkSkip = [img skipSizeForLoop:blkLp];

	proc = ^void(float *p, float *q, int skip, int len) {
		float	re, im;
		float	*bp, *bq;
        int		j, k, ix;

		bp = [img currentDataWithControl:blkLc];
		bq = bp + [img dataLength];
		for (k = 0; k < nBlk; k++) {
			re = im = 0;
			for (j = 0; j < blkLen; j++) {
				ix = (k * blkLen + j + ds) * skip;
				re += p[ix];
				im += q[ix];
			}
			bp[k * blkSkip] = re;
			bq[k * blkSkip] = im;
		}
    };
    [self applyComplex1dProc:proc forLoop:avgLp control:lc];

    [self copyIvarOf:img];

    return blkLp;
}

- (RecImage *)reorderWithDda:(int)dda phsLp:(RecLoop *)phsLp avgLp:(RecLoop *)avgLp phsArray:(int *)array
{
	RecImage		*img;
	RecLoopIndex	*phsLi, *avgLi;
	int				nPhs, nAvg;
	int				phs;
	int				srcLen, dstLen, siz;
	int				nRep = [[self topLoop] dataLength];
	RecLoopControl	*lc1, *lc2;		// src, dst
	float			*p1, *p2;
	float			*q1, *q2;
	int				i, j;
	int				o1, o2, d1, d2;

//	avgLp = [RecLoop loopWithDataLength:nAvg];
//	phsLp = [RecLoop loopWithDataLength:nPhs];
	nPhs = [phsLp dataLength];
	nAvg = [avgLp dataLength];
	img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:phsLp, avgLp, [self yLoop], [self xLoop], nil];

	srcLen = [self dataLength];
	dstLen = [img dataLength];

	lc1 = [self control];
	[lc1 deactivateXY];

	lc2 = [img control];
	phsLi = [lc2 loopIndexForLoop:phsLp];
	avgLi = [lc2 loopIndexForLoop:avgLp];
	siz = [img xDim] * [img yDim];

// protocol file format:
//	n_image lines, each having one delay number (float)
//	-100: dda, 0:off, positive numbers:delay in msec

	d1 = d2 = o1 = o2 = dda = 0;
	for (i = 0; i < nRep; i++) {
		phs = array[i];
		p1 = [self currentDataWithControl:lc1];	// src
		q1 = p1 + srcLen;
		if (phs == 0) {
			// stim off, add to phs 0 & 2
			if (i % 4 == 0) {
				// off 1 (phs 0)
				[phsLi setCurrent:0];
				[avgLi setCurrent:o1];
				o1++;	// index of avg loop
			} else {
				// off 2 (phs 2)
				[phsLi setCurrent:2];
				[avgLi setCurrent:o2];
				o2++;	// index of avg loop
			}
		} else
		if (phs == 1) {
			// delay 1 (phs 1)
			[phsLi setCurrent:1];
			[avgLi setCurrent:d1];
			d1++;	// index of avg loop
		} else {
			// delay 2  (phs 3)
			[phsLi setCurrent:3];
			[avgLi setCurrent:d2];
			d2++;	// index of avg loop
		}
		p2 = [img currentDataWithControl:lc2];
		q2 = p2 + [img dataLength];
		for (j = 0; j < siz; j++) {
			p2[j] = p1[j];
			q2[j] = q1[j];
		}
		if (d1 >= nAvg && d2 >= nAvg) break;
		[lc1 increment];
	}

	return img;
}

// ### not done yet
- (RecImage *)phsSliceAtIndex:(int)ix nPhs:(int)nphs dda:(int)dda
{
	RecImage	*img, *slice;
	RecLoop		*avgLp, *zLp;
	int			i, k, n, nImg;

	nImg = [self zDim];
	zLp = [self zLoop];
	n = 0;

	n = 0;
	for (i = 0; i < nImg; i++) {
		if ((i - dda) % nphs == ix) n++;	// ###
	}
	printf("n = %d\n", n);

	avgLp = [RecLoop loopWithDataLength:n];
	img = [RecImage imageOfType:[self type] withLoops:avgLp, [self yLoop], [self xLoop], nil];

	// copy slice
	k = 0;
	for (i = 0; i < nImg; i++) {
		if ((i - dda) % nphs == ix) {	// ###
			slice = [self sliceAtIndex:i forLoop:zLp];
			[img copySlice:slice atIndex:k];
			k++;
		}
	}

	return img;
}

// ### not done yet
- (RecImage *)phsSliceAtIndex:(int)ix withPhsArray:(int *)phsArray dda:(int)dda
{
	RecImage	*img, *slice;
	RecLoop		*avgLp, *zLp;
	int			i, k, n, nImg;

	nImg = [self zDim];
	zLp = [self zLoop];
	n = 0;
	for (i = 0; i < nImg; i++) {
		if (phsArray[i + dda] == ix) n++;	// ###
	}
	printf("n = %d\n", n);

	avgLp = [RecLoop loopWithDataLength:n];
	img = [RecImage imageOfType:[self type] withLoops:avgLp, [self yLoop], [self xLoop], nil];

	// copy slice
	k = 0;
	for (i = 0; i < nImg; i++) {
		if (phsArray[i + dda] == ix) {	// ###
			slice = [self sliceAtIndex:i forLoop:zLp];
			[img copySlice:slice atIndex:k];
			k++;
		}
	}

	return img;
}

// assuming [kz, phs, ky, kx]
- (void)driftCorr
{
    RecImage        *img, *ref, *corr, *param;

	img = [self avgForLoop:[self topLoop]];
//[img saveAsKOImage:@"IMG_drift"];
    ref = [img sliceAtIndex:0];
    corr = [img xyCorrelationWith:ref];
    param = [corr estShift];
    img = [self ftShiftBy:param];
//[img saveAsKOImage:@"IMG_drift_corr"];
    [self copyIvarOf:img];
}

- (RecImage *)reorderBoldWithDda:(int)dda nEvent:(int)ne evLen:(int)el slcLp:(RecLoop *)slcLp phsLp:(RecLoop *)phsLp avgLp:(RecLoop *)avgLp delay:(int)del
{
	RecImage		*img;
	RecLoopIndex	*slcLi, *phsLi, *avgLi;
	int				nPhs, nAvg;
	int				srcLen, dstLen, siz;
	RecLoopControl	*lc1, *lc2;		// src, dst
	float			*p1, *p2;		// src, dst
	int				i, j, k;
	int				slc, nSlc;
	int				ix_on, ix_off;

	img = [RecImage imageOfType:RECIMAGE_REAL withLoops:slcLp, phsLp, avgLp, [self yLoop], [self xLoop], nil];
	nSlc = [slcLp dataLength];
	nPhs = [phsLp dataLength];
	nAvg = [avgLp dataLength];
//[img dumpLoops];
	srcLen = [self dataLength];
	dstLen = [img dataLength];

	lc1 = [self control];
	[lc1 deactivateXY];
	lc2 = [img controlWithControl:lc1];
	slcLi = [lc2 loopIndexForLoop:slcLp];
	phsLi = [lc2 loopIndexForLoop:phsLp];
	avgLi = [lc2 loopIndexForLoop:avgLp];
	siz = [img xDim] * [img yDim];

	// event data
	for (slc = 0; slc < nSlc; slc++) {
		[lc1 rewind];
		ix_on = ix_off = 0;
		// skip dda
		for (i = 0; i < dda + del; i++) {
			[lc1 increment];
		}
		[slcLi setCurrent:slc];
		for (i = 0; i < ne; i++) {
			// off
			for (j = 0; j < el; j++) {
				p1 = [self currentDataWithControl:lc1];	// src
				[phsLi setCurrent:0];
				[avgLi setCurrent:ix_off];
				p2 = [img currentDataWithControl:lc2];
				for (k = 0; k < siz; k++) {
					p2[k] = p1[k];
				}
				ix_off++;
				[lc1 increment];
			}
			// on
			for (j = 0; j < el; j++) {
				p1 = [self currentDataWithControl:lc1];	// src
				[phsLi setCurrent:1];
				[avgLi setCurrent:ix_on];
				p2 = [img currentDataWithControl:lc2];
				for (k = 0; k < siz; k++) {
					p2[k] = p1[k];
				}
				ix_on++;
				[lc1 increment];
			}
		}
	}

	return img;
}

// select phs = ix1 & ix2 (time order is not changed)
- (RecImage *)makeSubSeriesWithDda:(int)dda phsLp:(RecLoop *)phsLp stim:(int)ix1 ref:(int)ix2
{
	RecImage		*img;
	RecLoop			*timeLp;
	int				i, j, len, siz, phs;
	int				nPhs;
	RecLoopControl	*lc;
	float			*p1, *p2;	// src, dst
	float			*q1, *q2;	// src, dst

	lc = [self control];
	[lc deactivateXY];
	len = [lc loopLength];
	len = Rec_down2Po2(len);
	siz = [self xDim] * [self yDim];
	nPhs = [phsLp dataLength];

	timeLp = [RecLoop loopWithDataLength:len * 2 / nPhs];

	img = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:timeLp, [self yLoop], [self xLoop], nil];
	p2 = [img data];
	q2 = p2 + [img dataLength];

	for (i = 0; i < dda * nPhs; i++) {
		[lc increment];
	}
	for (i = 0; i < len; i++) {
		phs = i % nPhs;
		if (phs == ix1 || phs == ix2) {
			p1 = [self currentDataWithControl:lc];
			q1 = p1 + [self dataLength];
			for (j = 0; j < siz; j++) {
				p2[j] = p1[j];
				q2[j] = q1[j];
			}
			p2 += siz;
			q2 += siz;
		}
		[lc increment];
	}
	return img;
}

// input img (self) is full data set including dda
- (void) addTestTag:(int)dda nPhs:(int)nPhs phsArray:(int *)phsArray amp:(float)amp	// put test marker
{
	int				i, len, kk;
	int				x, y;
	RecLoopControl	*lc;
	float			*p;

int dbg = 1;

	lc = [self control];
	[lc deactivateXY];
	len = [lc loopLength];

// phs/avg loop version
	if (dbg == 0) {
		for (i = 0; i < len; i++) {
			if (i >= dda) {
				p = [self currentDataWithControl:lc];
				p += [self dataLength];
				kk = (i - dda) % nPhs;
				if (kk == 1) {
					x = 35; y = 80;
					p[y * 128 + x + 1] += amp;
					p[y * 128 + x - 1] -= amp;
				} else
				if (kk == 3) {
					x = 40; y = 80;
					p[y * 128 + x + 1] += amp;
					p[y * 128 + x - 1] -= amp;
				}
			}
			[lc increment];
		}
	} else {
// phsArray version
		printf("add test tag using phsArray\n");
		int k1 = 0, k2 = 0;
		for (i = 0; i < len; i++) {	// len = nImages (530)
			p = [self currentDataWithControl:lc];
			p += [self dataLength];	// imag part
			kk = phsArray[i];
			if (kk == 1) {
				x = 35; y = 80;
				p[y * 128 + x + 1] += amp;
				p[y * 128 + x - 1] -= amp;
				k1++;
			} else
			if (kk == 2) {
				x = 40; y = 80;
				p[y * 128 + x + 1] += amp;
				p[y * 128 + x - 1] -= amp;
				k2++;
			}
			// else do nothing (dda)

			[lc increment];
		}
		printf("k1 = %d, k2 = %d\n", k1, k2);
	}
}

- (RecImage *)filtImag:(float)w forLoop:(RecLoop *)lp	// HPfilt imag (frac = 1)
{
	RecImage	*im;

	im = [self copy];
	[im takeImagPart];
	[im gauss1DHP:w forLoop:lp frac:1.0];

	return im;
}

- (void)singleFreqFilt:(int)freq forLoop:(RecLoop *)lp	// remove single frequency (single fourier coeff)
{
    void    (^proc)(float *p, float *q, int skip, int len);
	int		i, n;
	BOOL	cpx = NO;

	n = [lp dataLength];

    proc = ^void(float *p, float *q, int skip, int len) {
        int     i, ix;
		
        for (i = ix = 0; i < n; i++, ix += skip) {
			if (ix == n/2 + freq || ix == n/2 - freq
				|| ix == n/2 + freq * 3 || ix == n/2 - freq * 3) {
				p[ix] = 0;
				q[ix] = 0;
			}
        }
    };
	if ([self type] != RECIMAGE_COMPLEX) {
		cpx = NO;
		[self makeComplex];
	} else {
		cpx = YES;
	}
	[self fft1d:lp direction:REC_INVERSE];
	[self applyComplex1dProc:proc forLoop:lp];
	[self fft1d:lp direction:REC_FORWARD];
	if (!cpx) {
		[self takeRealPart];
	}
}


- (void)cycFilt:(int)len forLoop:(RecLoop *)lp	// remove single frequency (cyclic filter)
{
	float			*buf = (float *)malloc(sizeof(float) * len);
	float			*p;
	int				i, j, k;	// i, j: time, k:pixel
	int				skip;
	int				nOut, nPixel;
	RecLoopControl	*lc = [self control];

	[lc deactivateLoop:lp];
	nPixel = [lc loopLength];
	nOut = [lp dataLength] / len;
	skip = [self skipSizeForLoop:lp];

	for (k = 0; k < nPixel; k++) {
		// make template
		p = [self currentDataWithControl:lc];
		for (j = 0; j < len; j++) {
			buf[j] = 0;
		}
		for (i = 0; i < nOut; i++) {
			for (j = 0; j < len; j++) {
				buf[j] += p[(i * len + j) * skip];
			}
		}
		for (j = 0; j < len; j++) {
			buf[j] /= nOut;
		}
		// subtract template
		for (i = 0; i < nOut; i++) {
			for (j = 0; j < len; j++) {
				p[(i * len + j) * skip] -= buf[j];
			//	p[(i * len + j) * skip] = buf[j];
			}
		}
		[lc increment];
	}

	free(buf);
}

- (void)accumFilt:(int)len forLoop:(RecLoop *)lp	// accum filter
{
	float			*buf = (float *)malloc(sizeof(float) * len);
	float			*p, mn;
	int				i, j, k;	// i, j: time, k:pixel
	int				skip;
	int				nOut, nPixel;
	RecLoopControl	*lc = [self control];

	[lc deactivateLoop:lp];
	nPixel = [lc loopLength];
	nOut = [lp dataLength] / len;
	skip = [self skipSizeForLoop:lp];

	for (k = 0; k < nPixel; k++) {
		p = [self currentDataWithControl:lc];
		// accum
		for (i = 0; i < nOut; i++) {
			for (j = 0; j < len; j++) {
				buf[j] = p[(i * len + j) * skip];
			}
			mn = 0;
			for (j = 0; j < len; j++) {
				mn += buf[j];
			}
			mn /= len;
			for (j = 0; j < len; j++) {
	//			buf[j] -= mn;
			}
			buf[0] = 0;
			for (j = 1; j < len; j++) {
				buf[j] += buf[j - 1];
			}
			for (j = 0; j < len; j++) {
				p[(i * len + j) * skip] = buf[j];
			}
		}
		[lc increment];
	}
	free(buf);
}

- (void)leakyIntWithTau:(float)tau cutOff:(float)width forLoop:(RecLoop *)lp
{
	RecImage		*kern = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:lp, nil];
	float			*p;
	int				i, len = [lp dataLength];

// make kernel (exponential decay with time constant tau)
	p = [kern real];
	for (i = 0; i < len; i++) {
		p[i] = exp((float)i / tau);
	}
	[kern shift1d:lp];

// freq domain filter
	[self fft1d:lp direction:REC_INVERSE];
	[kern fft1d:lp direction:REC_INVERSE];
	[kern conjugate];
	[self multByImage:kern];
	if (width > 0) {
		[self fGauss1DHP:width forLoop:lp frac:1.0];
	}
	[self fft1d:lp direction:REC_FORWARD];
	[self takeRealPart];
}

- (void)cos1dForLoop:(RecLoop *)lp cyc:(int)cyc power:(int)pw lowPass:(BOOL)lpf
{
    [self fft1d:lp direction:REC_INVERSE];
    [self fCos1dForLoop:lp cyc:cyc power:pw lowPass:lpf];
    [self fft1d:lp direction:REC_FORWARD];
}

// cos 1d
- (void)fCos1dForLoop:(RecLoop *)lp cyc:(int)cyc power:(int)pw lowPass:(BOOL)lpf
{
    void    (^proc)(float *p, int n, int skip);
	float	*wt, th;
	int		i, n;

	n = [lp dataLength];
	wt = (float *)malloc(sizeof(float) * n);
	for (i = 0; i < n; i++) {
		th = ((float)i - n/2) * 2 * cyc * M_PI / n;
		if (lpf) {
			wt[i] = (cos(th) + 1) / 2;
		} else {
			wt[i] = (1 - cos(th)) / 2;
		}
		wt[i] = pow(wt[i], pw);
	}

    proc = ^void(float *p, int n, int skip) {
        int     i, ix;

        for (i = ix = 0; i < n; i++, ix += skip) {
            p[ix] *= wt[i];
        }
    };
    [self apply1dProc:proc forLoop:lp];
	free(wt);
}

- (RecImage *)removeDdaAndMakePo2:(int)dda
{
	int				i, n;
	RecLoopControl	*lc;
	RecImage		*img;
	float			*p1, *q1;
	float			*p2, *q2;

	lc = [self control];
	[lc deactivateXY];
	n = [lc loopLength];
	n = Rec_down2Po2(n);

	for (i = 0; i < dda; i++) {
		[lc increment];
	}

	img = [RecImage imageOfType:[self type] xDim:[self xDim] yDim:[self yDim] zDim:n];
	p1 = [self currentDataWithControl:lc];
	q1 = p1 + [self dataLength];
	p2 = [img data];
	q2 = p2 + [img dataLength];
	n *= [self xDim] * [self yDim];
	for (i = 0; i < n; i++) {
		p2[i] = p1[i];
		q2[i] = q1[i];
	}
	return img;
}

- (RecImage *)addPhs:(RecLoop *)phsLp protocol:(int)protocol dda:(int)dda	// reverse of above
{
	int				i, k, n, len;
	RecLoopControl	*lc;
	RecImage		*img;
	RecLoop			*avgLp;
	int				nAvg;
	float			*p, *q;

	lc = [self control];
	[lc deactivateXY];
	n = [lc loopLength];
	nAvg = (n - dda) / [phsLp dataLength];
	avgLp = [RecLoop loopWithDataLength:nAvg];
	img = [RecImage imageOfType:[self type] withLoops:avgLp, phsLp, [self yLoop], [self xLoop], nil];
	len = [self dataLength] - dda;
	p = [self data];
	q = [img data];

	for (k = 0; k < pixSize; k++) {
		for (i = 0; i < len; i++) {
			q[i] = p[i + dda];
		}
		p += [self dataLength];
		q += [img dataLength];
	}
	return img;
}

void
dump_coef(RecImage *coef)
{
	int		i, j;
	int		order, nLines;
	float	*p;

	order = [coef xDim];
	nLines = [coef outerLoopDim];
	p = [coef data];
	for (i = 0; i < nLines; i++) {
		printf("%d ", i);
		for (j = 0; j < order; j++) {
			printf("%f ", p[i * order + j]);
		}
		printf("\n");
	}
}

- (RecImage *)phsAvgForLoop:(RecLoop *)lp
{
	RecImage	*img;
	img = [self avgForLoop:lp];
	[img phase];
	return img;
}

- (RecImage *)phsVarForLoop:(RecLoop *)lp withMean:(RecImage *)mn
{
	RecImage	*img, *phs;

	img = [self copy];
	phs = [mn copy];
	[phs thToExp];
	[img cpxDivImage:phs];
	[img phase];
	[img square];
	img = [img avgForLoop:lp];

	return img;
}

- (RecImage *)phsSDForLoop:(RecLoop *)lp
{
	RecImage	*sd, *mn;

	mn = [self phsAvgForLoop:lp];
	sd = [self phsVarForLoop:lp withMean:mn];
	[sd sqroot];

	return sd;
}

- (RecImage *)sdMaskForLoop:(RecLoop *)lp andLoop:(RecLoop *)lp2 thres:(float)thres
{
	RecImage	*sd = [self phsSDForLoop:lp];

	sd = [sd avgForLoop:lp2];
	[sd thresAt:thres frac:NO];
	[sd logicalInv];

	return sd;
}

- (RecImage *)sdMaskForLoop:(RecLoop *)lp thres:(float)thres
{
	RecImage	*sd = [self phsSDForLoop:lp];
	[sd thresAt:thres frac:NO];
	[sd logicalInv];

	return sd;
}

- (void)subDataOf:(RecImage *)img
{
	int		i, len;
	float	*p, *q;

	len = [self dataLength] * [self pixSize];
	p = [self data];
	q = [img data];
	for (i = 0; i < len; i++) {
		p[i] -= q[i];
	}
}

// remove low freq phase (polynomial base)
- (void)baseline2Avg:(RecLoop *)avg phs:(RecLoop *)phs
{
	RecImage		*av, *msk;
	int				order = 5;
	RecLoop			*ordx, *ordy;
	RecImage		*coef;
	RecLoopControl	*lc;
	float			thres;
	
	thres = 0.30;	// 0.15

//	avg = [RecLoop findLoop:@"avg"];
//	phs = [RecLoop findLoop:@"phs"];
	av = [self avgForLoop:avg];
	av = [av avgForLoop:phs];
	msk = [self sdMaskForLoop:avg andLoop:phs thres:thres];
	[av maskWithImage:msk];
//[av saveAsKOImage:@"IMG_av_msk"];

	ordx = [RecLoop loopWithDataLength:order];
	ordy = [RecLoop loopWithDataLength:order];
	lc = [av control];
	[lc replaceLoop:[lc xLoop] withLoop:ordx];
	[lc replaceLoop:[lc yLoop] withLoop:ordy];
	coef = [RecImage imageOfType:RECIMAGE_REAL withControl:lc];
	[av pestPoly2d:coef];
//dump_coef(coef);

	[self pcorrPoly2d:coef];
}

// basically same as above... replace above when done
- (void)baseline3Rep:(RecLoop *)repLp
{
	RecImage		*av, *msk;
	int				order = 5;
	RecLoop			*ordx, *ordy;
	RecImage		*coef;
	RecLoopControl	*lc;
	float			thres;
	
	thres = 0.30;	// 0.15

	av = [self avgForLoop:repLp];
	msk = [self sdMaskForLoop:repLp thres:thres];
	[av maskWithImage:msk];
[av saveAsKOImage:@"IMG_av_msk"];

	ordx = [RecLoop loopWithDataLength:order];
	ordy = [RecLoop loopWithDataLength:order];
	lc = [av control];
	[lc replaceLoop:[lc xLoop] withLoop:ordx];
	[lc replaceLoop:[lc yLoop] withLoop:ordy];
	coef = [RecImage imageOfType:RECIMAGE_REAL withControl:lc];
	[av pestPoly2d:coef];
dump_coef(coef);

	[self pcorrPoly2d:coef];
}

// remove DC along zLoop (1d)
- (void)removeDC
{
	RecImage	*mn;

	mn = [self avgForLoop:[self zLoop]];
	[self subImage:mn];
}

// DSB (double side-band) decoding ### not correct yet
// can also interpreted as Gabor filter
- (RecImage *)dsbDecForLoop:(RecLoop *)lp width:(float)w center:(float)c
{
	RecImage	*imgp, *imgm;

	imgp = [self copy];
	[imgp fft1d:lp direction:REC_INVERSE];
	[imgp shift1d:lp by:(int)c];
	[imgp fGauss1DLP:w forLoop:lp];
	[imgp fft1d:lp direction:REC_FORWARD];
	imgm = [self copy];
	[imgm fft1d:lp direction:REC_INVERSE];
	[imgm shift1d:lp by:-(int)c];
	[imgm fGauss1DLP:w forLoop:lp];
	[imgm fft1d:lp direction:REC_FORWARD];
	[imgp addImage:imgm];

	return imgp;
}

// seems to be working, but slow
// 2D polynomial phase correction // ### not done yet (1d unwrap)
// ... coef changed to 2D (4/28) ... should be ok ...
// ### moved to nciRec.m (tmp)
- (void)pestPoly2d:(RecImage *)coef
{
	void			(^proc)(float *p, float *q, int xDim, int yDim);
	int				ordx, ordy;
	int				xDim, yDim;
	RecLoopControl	*lc_coef, *lc_self;

	ordx = [coef xDim];
	ordy = [coef yDim];
	xDim = [self xDim];
	yDim = [self yDim];

	lc_self = [self control];
	lc_coef = [coef controlWithControl:lc_self];

// image method version
    proc = ^void(float *p, float *q, int xDim, int yDim) {
		float *cp = [coef currentDataWithControl:lc_coef];
		Rec_est_poly_2d(cp, ordx, ordy, p, q, xDim, xDim);	// ### not done yet
    };
    [self applyComplex2dProc:proc control:lc_self];
}

// slow phase removal
- (void)baseline  // just reproduce nci_rec version of baseline()
{
    RecImage    *img;

//    img = [self gauss2d:3.0];
	img = [self copy];
	[img gauss2DLP:0.2];
    [img toUnitImage];
    [img conjugate];
    [self multByImage:img];
}

// quadrupole component removal
- (void)quadFilt:(float)qd
{
	RecImage	*wt;
	float		*p, x, y, r;
	int			i, j;
	int			xDim, yDim;

	[self fft2d:REC_INVERSE];

	xDim = [self xDim];
	yDim = [self yDim];
	wt = [RecImage imageOfType:RECIMAGE_REAL xDim:xDim yDim:yDim];
	p = [wt data];
	for (i = 0; i < yDim; i++) {
		y = ((float)i - yDim/2) / yDim;
		for (j = 0; j < xDim; j++) {
			x = ((float)j - xDim/2) / xDim;
			r = sqrt(x*x + y*y);
			p[i*xDim + j] = 1.0 / (r * qd + 1.0);
		}
	}
[wt saveAsKOImage:@"IMG_quad"];
	[self multByImage:wt];

	[self fft2d:REC_FORWARD];
}

// filter spike noise in phase image: thres in radian
- (void)spikeFilt:(float)thres
{
	RecLoopControl	*lc = [self control];
	int				i, j, n, loopLen;
	float			*p;
	
	[lc deactivateX];
	loopLen = [lc loopLength];
	[lc rewind];
	n = [self xDim];
	for (i = 0; i < loopLen; i++) {
		p = [self currentDataWithControl:lc];
		// ### 1d filter
		for (j = 0; j < n; j++) {
			if (fabs(p[j]) > thres) {
				p[j] = 0;
			}
		}
		
		[lc increment];
	}
}

void
dump_vector(NSString *path, float *p, int n)
{
	int		i;
	FILE	*fp = fopen([path UTF8String], "w");
	for (i = 0; i < n; i++) {
		fprintf(fp, "%d %g\n", i, p[i]);
	}
	fclose(fp);
}

// PCA type variance filter
// started on 3-7-2017 ###
// real (phase) only for the moment
// restarted on 3-14-2017 using SVD
// zLoop == timeLoop, real image
- (void)pcaFilt:(int)ns
{
	RecLoop			*timeLp;
	int				i, j, k, tDim, iDim;
	float			*p, *q, *qq;
//	NSString		*path;

	RecImage		*img, *eimg, *uimg, *err, *tmp;
	Num_mat			*A;
	Num_svd_result	*sres;
	Num_ica_result	*ires;
	int				nr, nc, n;
//	Num_mat			*W, *Sr;
	int				ncomp = 16;	// po2
	RecLoop			*ncLoop = [RecLoop loopWithDataLength:ncomp];

// input
	[self saveAsKOImage:@"IMG_pca_in"];

	timeLp = [self zLoop];
	tDim = [timeLp dataLength];
	iDim = [self xDim] * [self yDim];

// === SVD
	//make A matrix
	img = [RecImage imageOfType:RECIMAGE_REAL xDim:iDim yDim:tDim];
	q = [img data];
	p = [self data];
	for (i = 0; i < iDim * tDim; i++) {
		q[i] = p[i];
	}
	[img saveAsKOImage:@"IMG_A"];
	nr = [img yDim];
	nc = [img xDim];
	n = MIN(nr, nc);
	A = Num_im_to_m(img);
	// make U, s, Vt
//	U = Num_new_mat(nr, n);
//	Vt = Num_new_mat(n, nc);
//	s = Num_new_vec(n);
	// call LAPACK
TIMER_ST
	sres = Num_svd(A);
TIMER_END("svd")

	img = Num_v_to_im(sres->s);
	[img saveAsKOImage:@"IMG_s"];

	uimg = Num_m_to_im(sres->U);
	[uimg trans];
	[uimg replaceLoop:[uimg yLoop] withLoop:ncLoop offset:0];
	[uimg saveAsKOImage:@"IMG_U"];
	tmp = [uimg copy];
	[tmp fft1d:[tmp xLoop] direction:REC_INVERSE];
	[tmp saveAsKOImage:@"IMG_U_f"];
	
	img = Num_m_to_im(sres->Vt);
	[img saveAsKOImage:@"IMG_Vt"];

	// eigen image
//	eimg = [RecImage imageOfType:RECIMAGE_REAL xDim:128 yDim:128 zDim:tDim];	// -> preserve img dim
	eimg = [RecImage imageOfType:RECIMAGE_REAL xDim:[self xDim] yDim:[self yDim] zDim:tDim];
	[eimg copyImageData:img];
	[eimg saveAsKOImage:@"IMG_E"];
	tmp = [eimg kFilt:0];
	[tmp saveAsKOImage:@"IMG_E_dp"];

	tmp = [tmp oversample];
	tmp = [tmp oversample];
	[tmp saveAsKOImage:@"IMG_E_dp_hi"];
	
	Num_free_svd_result(sres);

//	ICA
	ires = Num_ica(A, ncomp);
	tmp = Num_m_to_im(ires->WK);
	eimg = [RecImage imageOfType:RECIMAGE_REAL xDim:[self xDim] yDim:[self yDim] zDim:ncomp];
	[eimg copyImageData:tmp];
	[eimg saveAsKOImage:@"IMG_map"];
	tmp = [eimg kFilt:0];
	[tmp saveAsKOImage:@"IMG_map_dp"];
	saveAsKOImage(ires->WX, @"IMG_out");
	tmp = Num_m_to_im(ires->WX);
	[tmp fft1d:[tmp xLoop] direction:REC_INVERSE];
	[tmp saveAsKOImage:@"IMG_out_f"];
	saveAsKOImage(ires->W, @"IMG_W");

	Num_free_ica_result(ires);

// ==== sort based on energy difference between noise and stim
// copy testProc / test4()


exit(0);	// 9-20-2017 ###

// make stim vector ===> caller should provide ###


// === filter based on SVD
// === update this part ### to consider stim vector component
// === calc correlation of noise vs stim
// === sort (optimize sorting scheme ###) 6-23-2017



	err = [RecImage imageOfType:RECIMAGE_REAL withImage:self];
	q = [uimg data];
	for (k = 0; k < ns; k++) {	// upto ns-th singular value
		// self -= E(i) x s(i)
		img = [eimg sliceAtIndex:k forLoop:[eimg zLoop]];
		p = [img data];
		qq = [err data];
		for (i = 0; i < tDim; i++) {
			for (j = 0; j < iDim; j++) {
				qq[j] = p[j] * q[i] * sres->s->data[k];
			}
			qq += iDim;
		}
		q  += tDim;
		// chk self & var
		[self subImage:err];
		printf("sval %d\n", k);
	}
	[err saveAsKOImage:@"IMG_err"];
	[self saveAsKOImage:@"IMG_flt"];

	Num_free_svd_result(sres);
}

// real only
- (void)outerProd:(RecImage *)wt andVector:(float *)bs len:(int)len	// ## make basis RecImage
{
	int		i, j;
	int		imgLen = [self xDim] * [self yDim];
	float	*p, *q;

	p = [wt data];
	q = [self data];
	for (i = 0; i < len; i++) {
		for (j = 0; j < imgLen; j++) {
			q[j] = p[j] * bs[i];
		}
		q += imgLen;
	}
}

// input: complex image (imag only) ===>  ##### real only complex !!! (imag = 0)
- (RecImage *)kFilt:(int)mmt
{
	[self makeComplex];

	switch (mmt) {
	default :
	case 0 :	// dipole
		return [self toDipole];
		break;
	case 1 :	// quadrupole
		return [self toQuadrupole];
		break;
	}
}

// moved to RecImage
/*
- (RecImage *)toDipole
{
    RecImage    *x, *y;

    x = [self copy];
    [x fft2d:REC_INVERSE];				// to freq
    y = [x copy];

    [x sinFilt:0 ang: 0.0];				// sin(th), odd func
    [x fft2d:REC_FORWARD];				// to space (pure imag)

    [y sinFilt:0 ang: M_PI/2.0];		// cos(th), odd func
    [y fft2d:REC_FORWARD];				// to space (pure imag)

    [x copyImagOf:y];					// make complex
    return x;
}
- (RecImage *)toQuadrupole
{
    RecImage    *x, *y;

    x = [self copy];
    [x fft2d:REC_INVERSE];				// to freq
    y = [x copy];

    [x sinFilt:1 ang: 0.0];				// sin(2*th), even func
    [x fft2d:REC_FORWARD];				// to space (pure real)

    [y sinFilt:1 ang: M_PI/4.0];		// cos(2*th), even func
    [y fft2d:REC_FORWARD];				// to space (pure real)

    [y copyRealOf:x];					// make complex
    return y;
}
*/

- (void)sinFilt:(int)mode ang:(float)ang
{
    void	(^proc)(float *p, float *q, int xDim, int yDim);

    proc = ^void(float *p, float *q, int xDim, int yDim) {
        // sin filt
		int		i, j, ix;
		float	x, y, th, w;

        for (i = ix = 0; i < yDim; i++) {
            y = (float)i - yDim/2;
            for (j = 0; j < xDim; j++, ix++) {
                x = (float)j - xDim/2;
				switch (mode) {
				default :
				case 0 :	// dipole
					th = atan2(y, x) - ang;  // was + ang
					break;
				case 1 :	// quadrupole
					th = atan2(y, x) * 2 - ang;  // was + ang
					break;
				}
				w = sin(th);
                p[ix] *= w;
                q[ix] *= w;
            }
        }
    };

    [self applyComplex2dProc:proc]; // main proc
}

// mg is real image
- (void)imagToPhaseUsingMag:(RecImage *)mg
{
    int     i;
    float   *p, *m;
    float   x;

    p = [self data];
    m = [mg data];
    for (i = 0; i < dataLength; i++) {
        if (m[i] == 0) {
            if (p[i] == 0) {
                x = 0.0;
            } else {
                x = 1.0;
            }
        } else {
            x = p[i] / m[i];
        }
        if (x > 1.0) x = 1.0;
        if (x < -1.0) x = -1.0;
        p[i] = asin(x);
    }
}

- (void)windowForLoop:(RecLoop *)lp mode:(int)filt width:(float)sd
{
    RecLoopControl  *lc;
    int             i, j, k, n, ix;
    int             dim, skip;
    float           *p;
    float           *w, x;

	lc = [self control];
	[lc deactivateLoop:lp];
	n = [lc loopLength];
    dim = [lp dataLength];
    skip = [self skipSizeForLoop:lp];

    w = (float *)malloc(dim * sizeof(float));
    for (j = 0; j < dim; j++) {
        x = (float)(j - dim/2) * 2 / dim;
        switch (filt) {
        case -1 :
            w[j] = 1.0;
            break;
        case 0 :    // cos
            w[j] = cos(x * 0.5 * M_PI / sd);
            break;
        case 1 :    // tri
            if (x > 0) {
                w[j] = 1 - x;
            } else {
                w[j] = 1 + x;
            }
            break;
        case 2 :    // gauss
            w[j] = exp(-3 * x * x / (sd*sd));
            break;
        }
    }

    for (i = 0; i < n; i++) {
        p = [self currentDataWithControl:lc];
        for (k = 0; k < pixSize; k++) {
            for (j = ix = 0; j < dim; j++, ix+=skip) {
                p[ix] *= w[j];
            }
            p += dataLength;
        }
        [lc increment];
    }

    free(w);
}

//	####
- (RecImage *)meanDiffWithRef:(RecImage *)ref forLoop:(RecLoop *)lp
{
	RecImage	*m_ref,  *m_tgt;

	m_ref = [ref avgForLoop:lp];
	m_tgt = [self avgForLoop:lp];
	[m_tgt subImage:m_ref];
	return m_tgt;
}

- (RecImage *)sdDiffWithRef:(RecImage *)ref forLoop:(RecLoop *)lp
{
	RecImage	*m_ref,  *m_tgt;
	RecImage	*sd_ref, *sd_tgt;

	m_ref = [ref avgForLoop:lp];
	sd_ref = [ref sdForLoop:lp withMean:m_ref];
//[sd_ref saveAsKOImage:@"IMG_sdref"];
	m_tgt = [self avgForLoop:lp];
	sd_tgt = [self sdForLoop:lp withMean:m_tgt];
//[sd_tgt saveAsKOImage:@"IMG_sdtgt"];
	[sd_tgt subImage:sd_ref];
//[sd_tgt saveAsKOImage:@"IMG_sddif"];
	return sd_tgt;
}

- (RecImage *)difForLoop:(RecLoop *)phs stim:(int)ix1 ref:(int)ix2
{
	RecImage	*df, *stim, *ref;

	stim = [self sliceAtIndex:ix1 forLoop:phs];
	ref  = [self sliceAtIndex:ix2 forLoop:phs];

	df = [stim copy];
	[df subImage:ref];

	return df;
}

- (RecImage *)difForLoop:(RecLoop *)phs stim:(int)ix1 ref1:(int)ix2 ref2:(int)ix3
{
	RecImage	*df, *stim, *ref1, *ref2;

	stim = [self sliceAtIndex:ix1 forLoop:phs];
	ref1  = [self sliceAtIndex:ix2 forLoop:phs];
	ref2  = [self sliceAtIndex:ix3 forLoop:phs];
	[ref1 addImage:ref2];
	[ref1 multByConst:0.5];

	df = [stim copy];
	[df subImage:ref1];

	return df;
}

// "differentiate" along loop (self is not altered)
// ### not working yet ????###
- (RecImage *)dxForLoop:(RecLoop *)lp
{
	RecImage		*res;
	RecLoopControl	*lc;
	int				i, j, ix, n, skip, loopLen;
	float			*p, *pp;

	res = [RecImage imageWithImage:self];
	n = [lp dataLength];
	loopLen = [self dataLength] / n;
	skip = [self skipSizeForLoop:lp];
	lc = [res control];
	[lc deactivateLoop:lp];
	for (i = 0; i < loopLen; i++) {
		p = [self currentDataWithControl:lc];
		pp = [res currentDataWithControl:lc];
		for (j = ix = 0; j < n; j++, ix += skip) {
			if (j == 0) {
				pp[ix] = 0;
			} else {
				pp[ix] = p[ix] - p[ix - skip];
			}
		}
	// complex ###
		[lc increment];
	}
	return res;
}

// not done yet ###
- (void)varFiltForLoop:(RecLoop *)lp thres:(float)th
{
	RecImage		*mg, *sd;
	RecLoopControl	*srcLc, *dstLc;
	int				i, j, n, loopLen;
	float			*p, *q, *pp, *qq;

	mg = [self copy];
	[mg magnitude];

	sd = [self sdForLoop:lp];	// central sd

	[sd multByConst:th];

	dstLc = [mg control];
	[dstLc deactivateXY];
	loopLen = [dstLc loopLength];
	srcLc = [RecLoopControl controlWithControl:dstLc forImage:sd];
	n = [mg xDim] * [mg yDim];
	for (i = 0; i < loopLen; i++) {
		pp = [self currentDataWithControl:dstLc];
		qq = pp + [self dataLength];
		p = [mg currentDataWithControl:dstLc];
		q = [sd currentDataWithControl:srcLc];
	//	printf("p = %0lx\n", (long)p);
	//	printf("q = %0lx\n", (long)q);
		for (j = 0; j < n; j++) {
			if (p[j] > q[j]) {
				pp[j] = qq[j] = 0;
			}
		}
		[dstLc increment];
	}
}

// F-stat not correct yet ###
// t-stat not implemented yet ###)
- (RecImage *)pImageWithRef:(RecImage *)ref complex:(BOOL)cpx
{
	if (cpx) {
		return [self pImageCpxWithRef:(RecImage *)ref];
	} else {
		// ### t-stat not implemented yet
		return nil;
	}
}

// test of mean diff of cpx numbers
- (RecImage *)pImageCpxWithRef:(RecImage *)ref
{
	RecImage	*pimg;
	RecImage	*mn1, *mn2;
	RecImage	*var1, *var2;
	RecImage	*ser;
	int			n1, n2;

	n1 = [self zDim];
	n2 = [ref  zDim];
	printf("n1 / n2 = %d / %d\n", n1, n2);
	mn1 = [self avgForLoop:[self zLoop]];
	mn2 = [ref  avgForLoop:[ref  zLoop]];

	var1 = [self varForLoop:[self zLoop] withMean:mn1];
	var2 = [ref  varForLoop:[ref  zLoop] withMean:mn2];

	[var1 multByConst:n1];
	[var2 multByConst:n2];

//[var1 sqroot];
//[var2 sqroot];
//dump_hist(var1, 100, "hist1.dat");
//dump_hist(var2, 100, "hist2.dat");
//exit(0);


	[var1 addImage:var2];
	[var1 multByConst:1.0/(n1 + n2 - 2.0)];

	ser = [var1 multByConst:(1.0/n1 + 1.0/n2)];	// ff
//	ser = [var1 multByConst:(1.0/(n1 + n2))];	// ff


//[ser saveAsKOImage:@"IMG_ser"];
	[mn1 subImage:mn2];
	[mn1 square];	// complex

	pimg = [mn1 copy];
	[pimg divImage:ser];	// F value

//dump_tdist(1);

dump_hist(pimg, 100, "hist1.dat");
exit(0);




//[pimg saveAsKOImage:@"IMG_fpval"];
	[pimg f2pDf1:2 df2:(n1 + n2 - 2)];	// 2
//	[pimg f2pDf1:n1 df2:n2];

	return pimg;
}

// assumes [slc, phs, avg, y, x]
- (RecImage *)pImageWithStim:(int)ix1 ref:(int)ix0 complex:(BOOL)cpx	// t or F-stat
{
	RecImage	*s1, *s0;
	RecImage	*m1, *m0;	// 1:stim, 0:ref
	RecImage	*v1, *v0;	// 1:stim, 0:ref
	RecLoop		*avg, *phs;
	int			n1, n0, df;		// for this method, n1 == n0 == [avg dataLength]

	phs = [[self loops] objectAtIndex:1];
	avg = [[self loops] objectAtIndex:2];
	n1 = n0 = [avg dataLength];
	
	s1 = [self sliceAtIndex:ix1 forLoop:phs];
	s0 = [self sliceAtIndex:ix0 forLoop:phs];

	m1 = [s1 avgForLoop:avg];
	m0 = [s0 avgForLoop:avg];
	v1 = [s1 varForLoop:avg withMean:m1];
	v0 = [s0 varForLoop:avg withMean:m0];
	df = n1 + n0 - 2;

	[m1 subImage:m0];

//[m1 saveAsKOImage:@"IMG_m1"];
//[v1 saveAsKOImage:@"IMG_v1"];
//[v0 saveAsKOImage:@"IMG_v0"];

	[v1 multByConst:n1 - 1];
	[v0 multByConst:n0 - 1];
	[v1 addImage:v0];
	[v1 multByConst:1.0/df];
	[v1 multByConst:(1.0/n1 + 1.0/n0)];
	[v1 sqroot];
	[m1 divImage:v1];
	
//[v1 saveAsKOImage:@"IMG_ser"];
//[m1 saveAsKOImage:@"IMG_timg"];

	[m1 t2pDf:df];
//	[m1 pvalThres:0.01];
//[m1 saveAsKOImage:@"IMG_pimg"];

	return m1;
}

- (void)f2pDf1:(int)df1 df2:(int)df2
{
	int		i, n = [self dataLength];
	float	*p = [self data];
	float	pval;

	for (i = 0; i < n; i++) {
		pval = Num_f2p(p[i], df1, df2);
		if (pval > 0) {
			pval = -log10(pval);
		} else {
			pval = 0;
		}
		p[i] = pval;
	}
}

- (void)t2pDf:(int)df
{
	int		i, n = [self dataLength];
	float	*p = [self data];
	float	pval;

	for (i = 0; i < n; i++) {
		pval = Num_t2p(p[i], df);
		if (pval > 0) {
			pval = -log10(pval);
		} else {
			pval = 0;
		}
		p[i] = pval;
	}
}

- (void)pvalThres:(float)th
{
	int		i, n = [self dataLength];
	float	*p = [self data];
	float	thres = -log(th);

	for (i = 0; i < n; i++) {
		if (p[i] < thres) {
			p[i] = 0;
		}
	}
}

// === NumKit version ===
/*
void
Num_ctmean(float *p1, float *q1, int n1, float *p2, float *q2, int n2, float *t2val, float *pval, float *th)
{
    float   var, serr2, t2, mr1, mi1, mr2, mi2, v1, v2;
    int     df;

    Num_cavevar(p1, q1, n1, &mr1, &mi1, &v1);
    Num_cavevar(p2, q2, n2, &mr2, &mi2, &v2);

    df = n1 + n2 - 2;
    var = (v1 * (n1 - 1) + v2 * (n2 - 1)) / df;
    serr2 = var * (1.0 / n1 + 1.0 / n2);

    mr1 -= mr2; mi1 -= mi2;
    t2 = mr1 * mr1 + mi1 * mi1;
    t2 /= serr2;

    *t2val = t2;
    *pval = Num_f2p(t2, 2, 2 * n1 - 1) * 2;
 //   *pval = Num_f2p(t2, 2, n1 - 1);
    *th = atan2(mi1, mr1);
}
*/

// returns complex mean / var test p-value
// mag : -log10(p), phs : phase of dipole
- (RecImage *)pImageWithRef:(RecImage *)ref forLoop:(RecLoop *)avgLp thres:(float)thres mode:(int)mode
{
    RecImage        *img;
    RecLoopControl  *srcLc, *dstLc, *refLc, *outerLc;
    RecLoop         *phsLp, *zLp;
    float           *p, *q;
    float           *pp, *qq;
    float           *re1, *re2;
    float           *im1, *im2;
    float           dr, di, th = 0;
    float           tval = 0, pval = 0;
    int             i, loopLen, dataLen;
    int             j, n, srcSkip, refSkip;
    int             k, outerLen;

    n = [avgLp dataLength];
    if (n < 4) {
        printf("navg too small\n");
        return nil;
    }

    re1 = VECTOR(n); //p1->re;
    im1 = VECTOR(n); //p1->im;
    re2 = VECTOR(n); //p2->re;
    im2 = VECTOR(n); //p2->im;

    srcLc = [RecLoopControl controlForImage:self];  // [phs, avg, z, y, x]
    refLc = [RecLoopControl controlWithControl:srcLc forImage:ref];
    dstLc = [RecLoopControl controlWithControl:srcLc];
    [dstLc removeLoop:avgLp];                       // [phs,      z, y, x]
    img = [RecImage imageOfType:RECIMAGE_COMPLEX withControl:dstLc];

    loopLen = [self xDim] * [self yDim] * [self zDim];
    dataLen = [img dataLength];
    outerLc = [RecLoopControl controlWithControl:dstLc];
    [outerLc deactivateXYZ];
    [outerLc deactivateLoop:avgLp];
    outerLen = [outerLc loopLength];
    [srcLc rewind];
    [srcLc deactivateLoop:avgLp];
    srcSkip = [self skipSizeForLoop:avgLp];
    refSkip = [ref skipSizeForLoop:avgLp];      // usually same as above, but depends on loop order

    for (i = 0; i < loopLen; i++) {
        for (k = 0; k < outerLen; k++) {
			switch (mode) {
			case NCI_REAL :
				p = [self currentDataWithControl:srcLc];
				for (j = 0; j < n; j++) {
					re1[j] = *p;
					p += srcSkip;
				}
				// copy 1d ref
				p = [ref currentDataWithControl:refLc];
				for (j = 0; j < n; j++) {
					re2[j] = *p;
					p += refSkip;
				}
				Num_ttest(re1, n, re2, n, &tval, &pval);
				break;
			case NCI_CPX_MEAN :
				// copy 1d data
				p = [self currentDataWithControl:srcLc];
				q = p + [self dataLength];
				for (j = 0; j < n; j++) {
					re1[j] = *p;
					im1[j] = *q;
					p += srcSkip;
					q += srcSkip;
				}
				// copy 1d ref
				p = [ref currentDataWithControl:refLc];
				q = p + [ref dataLength];
				for (j = 0; j < n; j++) {
					re2[j] = *p;
					im2[j] = *q;
					p += refSkip;
					q += refSkip;
				}
				Num_ctmean(re1, im1, n, re2, im2, n, &tval, &pval, &th);
				break;
			case NCI_CPX_VAR :
				// copy 1d data
				p = [self currentDataWithControl:srcLc];
				q = p + [self dataLength];
				for (j = 0; j < n; j++) {
					re1[j] = *p;
					im1[j] = *q;
					p += srcSkip;
					q += srcSkip;
				}
				// copy 1d ref
				p = [ref currentDataWithControl:refLc];
				q = p + [ref dataLength];
				for (j = 0; j < n; j++) {
					re2[j] = *p;
					im2[j] = *q;
					p += refSkip;
					q += refSkip;
				}
				Num_ctvar(re1, im1, n, re2, im2, n, &tval, &pval, &th);
				break;
			}
			if (pval > 0) {
				pval = -log10(pval);
			} else {
				pval = 0;
			}
			if (mode == NCI_REAL) {
				pp = [img currentDataWithControl:dstLc];
				*pp = pval;
			} else {
				if (pval < -log10(thres)) pval = 0;
				dr = pval * cos(th);
				di = pval * sin(th);
				// set pixel val
				pp = [img currentDataWithControl:dstLc];
				qq = pp + dataLen;
				*pp = dr;
				*qq = di;
			}
            [outerLc increment];
        }
        [srcLc increment];  // XYZ
    }
    // swap phs <-> z
    zLp = [img zLoop];
    phsLp = [img topLoop];
    if (![zLp isEqual:phsLp]) {
        [img swapLoop:phsLp withLoop:zLp];
    }

	free(re1); free(im1);
	free(re2); free(im2);

    return img;
}

// mode: 0x00:sum, 0x10:over, 0x00:p-mag, 0x01:complex
- (RecImage *)fusionWith:(RecImage *)sg gain:(float)gain mode:(int)mode
{
    RecImage        *img;
    RecLoopControl *lc, *selfLc, *sgLc;
    int             i, imgLen, sgLen;
    float           *p, *q, *pp;
    float           *r, *g, *b;
    float           rho, th;
    float           rr, gg, bb;
    float           ofs = 0 * M_PI;


    img = [RecImage imageOfType:RECIMAGE_COLOR withImage:sg];
    lc = [RecLoopControl controlForImage:img];
    selfLc = [RecLoopControl controlWithControl:lc forImage:self];
	sgLc = [RecLoopControl controlWithControl:lc forImage:sg];

    sgLen = [sg dataLength];
	imgLen = [img dataLength];

    [lc rewind];
    for (i = 0; i < imgLen; i++) {
        pp = [self currentDataWithControl:selfLc];
        p = [sg currentDataWithControl:sgLc];
        q = p + sgLen;

        r = [img currentDataWithControl:lc];
        g = r + imgLen;
        b = g + imgLen;

        rho = p[0] * p[0] + q[0] * q[0];
        rho = sqrt(rho);
        th = atan2(q[0], p[0]) + ofs;

        switch (mode & 0x0f) {
        case 0 : // p-mag        //    if (rho > 1.0) {            // p < 0.01
            rr = gain * rho;
            gg = gain * rho;
            bb = 0;
            break;
        case 1 : // p-dir
        default :
			// make common func (for ref col etc) ###
        //    rr = gain * rho * (0.2 + cos(th));
        //    gg = gain * rho * (0.1 + cos(th + M_PI * 2 / 3));
        //    bb = gain * rho * (0.4 + cos(th + M_PI * 4 / 3));
            rr = gain * rho * (1 + cos(th));
            gg = gain * rho * (1 + cos(th + M_PI * 2 / 3));
            bb = gain * rho * (1 + cos(th + M_PI * 4 / 3));
            break;
        }

        switch (mode & 0x10) {
        case 0x00 :    // sum
            r[0] = pp[0] + rr;
            g[0] = pp[0] + gg;
            b[0] = pp[0] + bb;
            break;
        case 0x10 :    // over
        default :
            if (rho > 0) {
                r[0] = rr;
                g[0] = gg;
                b[0] = bb;
            } else {
                r[0] = g[0] = b[0] = pp[0];
            }
            break;
        }
        [lc increment];
    }

    return img;
}

// self is mag
- (RecImage *)fusionWithPimg:(RecImage *)sg gain:(float)gain
{
    RecImage        *img;
    int             i, len, ix, lutLen = 100;
	int				st, wd;
    float           *p, *pp;
    float           *r, *g, *b;
	float			*rtab, *gtab, *btab;

    img = [RecImage imageOfType:RECIMAGE_COLOR withImage:self];
	len = [img dataLength];
	rtab = (float *)malloc(sizeof(float) * lutLen);
	gtab = (float *)malloc(sizeof(float) * lutLen);
	btab = (float *)malloc(sizeof(float) * lutLen);

	pp = [self data];
	p  = [sg data];

	r = [img data];
	g = r + len;
	b = g + len;

// lut (make gen function)
	st = 0;
	wd = 20;	// 0.2 n
	for (i = 0; i < wd; i++) {
		ix = i + st;
		rtab[ix] = sin(0.5 * i * M_PI / wd);
		gtab[ix] = 0;
		btab[ix] = 0;
	}
	st = 20;
	wd = 60;	// 0.6 n
	for (i = 0; i < wd; i++) {
		ix = i + st;
		rtab[ix] = 1;
		gtab[ix] = sin(0.5 * i * M_PI / wd);
		btab[ix] = 0;
	}
	st = 80;
	wd = 20;	// 0.2 n
	for (i = 0; i < wd; i++) {
		ix = i + st;
		rtab[ix] = 1;
		gtab[ix] = 1;
		btab[ix] = sin(0.5 * i * M_PI / wd);
	}

// fusion
    for (i = 0; i < len; i++) {
		ix = (int)p[i] * 100 / 20.0;
		if (ix < 0) {
			ix = 0;
		}
		if (ix >= lutLen) {
			ix = lutLen - 1;
		}
		r[i] = pp[i] + (rtab[ix] - 1.0) * gain;
		g[i] = pp[i] + (gtab[ix] - 1.0) * gain;
		b[i] = pp[i] + (btab[ix] - 1.0) * gain;
    }

	free(rtab);
	free(gtab);
	free(btab);
    return img;
}

// q & d version ...
- (void)clearSliceAtIndex:(int)ix
{
	int				i, n = [self xDim] * [self yDim];
	RecLoop			*tLp = [self zLoop];
	RecLoopControl	*lc = [self control];
	RecLoopIndex	*li = [lc loopIndexForLoop:tLp];
	float			*p, *q;

	[li setCurrent:ix];
	p = [self currentDataWithControl:lc];
	q = p + [self dataLength];
	for (i = 0; i < n; i++) {
		p[i] = 0;
		q[i] = 0;
	}
}

/*
- (RecImage *)cvarForLoop:(RecLoop *)lp withMean:(RecImage *)mn
{
    RecImage            *img;
    void                (^proc)(float *dp, float *dq, float *sp, float *sq, int len, int skip);
    int                 n;
 
	[self subImage:mn];

    n = [lp dataLength];
    proc = ^void(float *dp, float *dq, float *sp, float *sq, int len, int skip) {
        int     j, ix;
        float   v;
		float	re, im;

// variance
		v = 0;
		for (j = ix = 0; j < len; j++, ix += skip) {
			re = sp[ix];
			im = sq[ix];
			v += re * re + im * im;
		}
        *dp = v / (len - 1);;
        *dq = 0;
    };

	img = [self applyComplexProjProc:proc forLoop:lp];

    return img;
}

- (RecImage *)cvarForLoop:(RecLoop *)lp
{
    RecImage            *img;
    void                (^proc)(float *dp, float *dq, float *sp, float *sq, int len, int skip);
    int                 n;
 
    n = [lp dataLength];
    proc = ^void(float *dp, float *dq, float *sp, float *sq, int len, int skip) {
        int     j, ix;
        float   mr, mi, v;
		float	re, im;

		vDSP_meanv(sp, skip, &mr, len);
		vDSP_meanv(sq, skip, &mi, len);

// variance
		v = 0;
		for (j = ix = 0; j < len; j++, ix += skip) {
			re = sp[ix]     - mr;
			im = sq[ix] - mi;
			v += re * re + im * im;
		}
        *dp = v / (len - 1);;
        *dq = 0;
    };

	img = [self applyComplexProjProc:proc forLoop:lp];

    return img;
}
*/

// color ref
- (RecImage *)colMapImage:(int)mode
{
	RecImage	*img;
	int			i, j, xDim, yDim, len, ix;
	float		*rr, *gg, *bb;
	float		x, y;

	xDim = [self xDim];
	yDim = [self yDim];
	img = [RecImage imageOfType:RECIMAGE_COLOR xDim:xDim yDim:yDim];
	len = [img dataLength];
	rr = [img data];
	gg = rr + len;
	bb = gg + len;
	for (i = 0; i < yDim; i++) {
		y = ((float)i - yDim/2) / yDim;
		for (j = 0; j < xDim; j++) {
			x = ((float)j - xDim/2) / xDim;
//			r = sqrt(x*x + y*y);
//			th = atan2(y, x);
			ix = i * xDim + j;
			// make common func (for ref col etc) ###
			nciCalcColor(mode, x, y, rr + ix, gg + ix, bb + ix);
		}
	}
	return img;
}

// remove const offset from each slice
- (void)removeDC2dWithMask:(RecImage *)mask
{
	RecImage	*slice;
	int			i, j, n, nSlice, npix;
	float		mn, *p, *mk;

	nSlice = [self zDim];
	mk = [mask data];
	for (i = 0; i < nSlice; i++) {
		slice = [self sliceAtIndex:i forLoop:[self zLoop]];
		n = [slice dataLength];
		p = [slice data];
		mn = 0;
		npix = 0;
		for (j = 0; j < n; j++) {
			if (mk[j] > 0) {
				mn += p[j];
				npix += 1;
			}
		}
		mn /= npix;
		for (j = 0; j < n; j++) {
			if (mk[j] > 0) {
				p[j] -= mn;
				p[j] *= mk[j];
			}
		}
		[self copySlice:slice atIndex:i];
	}
}

// #### not done yet ... w not working ?
- (void)gammaFiltForLoop:(RecLoop *)lp width:(float)w
{
	int			i, nf = [lp dataLength];
	int			nt = 64, fct = 4;
	RecImage	*kern, *fkern;
	float		*p, *q, x;

	nt *= fct;
	w *= fct;
	kern = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nt];
	p = [kern real];
	for (i = 0; i < nt; i++) {
		x = (float)i;
		p[i] = (2 * x - x * x / w) * exp(-x / w);
	}

	nf *= fct;
	fkern = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:nf];
	q = [fkern real];
	for (i = 0; i < nt; i++) {
		q[i] = p[i];
	}
	[fkern shift1d:[fkern xLoop]];
[fkern saveAsKOImage:@"IMG_1"];
	
	[fkern fft1d:[fkern xLoop] direction:REC_INVERSE];
	[fkern crop:[fkern xLoop] to:[fkern xDim]/fct];
	[fkern conjugate];
[fkern saveAsKOImage:@"IMG_3"];

	[self fft1d:lp direction:REC_INVERSE];
	[self multByImage:fkern];
	[self fft1d:lp direction:REC_FORWARD];
}

// BOLD response curve
- (RecImage *)genRespWithLoop:(RecLoop *)lp dt:(float)dt
{
	RecImage	*img;
	int			i, n;
	float		x, a;
	float		*p;
	// model
	float		pd = 1.0;
	float		d1 = 2.5;	// sec
	float		d2 = 3.5;	// sec
	float		t1 = 5.0;	// sec
	float		t2 = 2.0;	// sec

	img = [RecImage imageOfType:RECIMAGE_REAL withLoops:lp, nil];
	p = [img data];
	n = [lp dataLength];
	for (i = 0; i < n; i++) {	// x 10msec
		x = dt * i;	// sec
		if (x >= 0 && x < d1) {
			a = 0;
		} else
		if (x < 10 + d2) {
			a = 1.0 - exp(-(x - d1) / t1);
		} else {
			a = 1.0 - exp(-(10.0 + d2 - d1) / t1);
			a *= exp(-(x - 10 - d2) / t2);
		}
		a = pd * (a - 0.4);
		p[i] = a;
	//	printf("%f %f\n", x, a);
	}
	
	return img;
}

// debug
- (void)dumpPhaseSDForLoop:(RecLoop *)lp image:(NSString *)path histo:(NSString *)histo
{
	float		*x, *p;
	int			i, n = 500;
	RecImage	*sd = [self phsSDForLoop:lp];
	FILE		*fp;

	[sd saveAsKOImage:path];
//	n = [lp dataLength];

	if (histo) {
		p = (float *)malloc(sizeof(float) * n);
		x = (float *)malloc(sizeof(float) * n);
		[sd histogram:p x:x min:0 max:2.0 binSize:n filt:NO];

		fp = fopen([histo UTF8String], "w");
		for (i = 0; i < n; i++) {
			fprintf(fp, "%f %f\n", x[i], p[i]);
		}
		free(p);
		free(x);
		fclose(fp);
	}
}

@end

void
nciCalcColor(int mode, float x, float y, float *r, float *g, float *b)
{
	float		rd, th;

	switch (mode) {
	case 0 :	// direction
		rd = sqrt(x*x + y*y);
		th = atan2(y, x);
		if (rd< 0.2 || rd > 0.5) {
			*r = 0.5;
			*g = 0.5;
			*b = 0.5;
		} else {
			*r = cos(th);
			*g = cos(th + M_PI * 2 / 3);
			*b = cos(th + M_PI * 4 / 3);
		}
		break;
	case 1 :	// direction x 2
		rd = sqrt(x*x + y*y);
		th = atan2(y, x);
		if (rd< 0.2 || rd > 0.5) {
			*r = 0.5;
			*g = 0.5;
			*b = 0.5;
		} else {
			*r = cos(th * 2);
			*g = cos(th * 2 + M_PI * 4 / 3);
			*b = cos(th * 2 + M_PI * 8 / 3);
		}
		break;
	}
}

// 2D phase est (tmp, moved from RecKit/RecUtil.m)
void
Rec_est_poly_2d(float *coef, int ordx, int ordy, float *re, float *im, int xDim, int yDim)
{
//	RecChebSetup	*setup;
	RecImage		*kern;
	Num_param		*param;
	float			*mg, *ph, *est;
	int				i, len, pdim;
	float			(^cost)(float *);
	int				iter;
	float			err;

	kern = Rec_cheb_mk_2d_kern(xDim, yDim, ordx, ordy);

	len = xDim * yDim;
	pdim = ordx * ordy;

	mg = (float *)malloc(sizeof(float) * len);	// magnitude & mask
	ph = (float *)malloc(sizeof(float) * len);
	est = (float *)malloc(sizeof(float) * len);

    cost = ^float(float *prm) {
		float	cst;
		int		ii;
		// calc ms error over image (within mask)
		Rec_cheb_2d_expansion(est, prm, kern);
		for (ii = 0; ii < len; ii++) {
			if (mg[ii] == 0) {
				est[ii] = ph[ii];	// dif = 0
			}
		}
		cst = Rec_L2_dist(ph, est, len);
        return cst;    
    };

	// calc mg and phs
	for (i = 0; i < len; i++) {
		mg[i] = sqrt(re[i] * re[i] + im[i] * im[i]);
		if (re[i] == 0) {
			mg[i] = ph[i] = 0;
		} else {
			ph[i] = atan2(im[i], re[i]);
		}
	}
	// 1d unwrap
	if (1) {
		for (i = 0; i < yDim; i++) {
			Rec_unwrap_1d(ph + (i * xDim), xDim, 1);
		}
		for (i = 0; i < xDim; i++) {
			Rec_unwrap_1d(ph + i, xDim, xDim);
		}
	}

// chk input
/*
{
	RecImage	*tmp;
	float		*p;
	tmp = [RecImage imageOfType:RECIMAGE_REAL xDim:xDim yDim:yDim];
	p = [tmp data];
	for (i = 0; i < len; i++) {
		p[i] = ph[i];
	}
	[tmp saveAsKOImage:@"IMG_ph_unwrap"];
}
*/
	// powell
	param = Num_alloc_param(pdim);
	for (i = 0; i < pdim; i++) {
		param->data[i] = 0;
	}

//TIMER_ST
	iter = Num_powell(param, cost, &err);
//TIMER_END("powel")
printf("2d pest: iter = %d, err = %f\n", iter, err);
	// result
	for (i = 0; i < pdim; i++) {
		coef[i] = param->data[i];
	}

	// cleanup
	Num_free_param(param);
	free(mg);
	free(ph);
	free(est);
}

void
dump_hist(RecImage *img, int n, char *path)
{
	FILE	*fp;
	float	*hist, *x, mx, mn;
	int		i;

	hist = (float *)malloc(sizeof(float) * n);
	x    = (float *)malloc(sizeof(float) * n);
	mx = [img maxVal];
	mn = [img minVal];
	[img histogram:hist x:x min:mn max:mx binSize:n filt:NO];
	fp = fopen(path, "w");
	for (i = 0; i < n; i++) {
		fprintf(fp, "%f %f\n", x[i], hist[i]);
	}
	fclose(fp);
	free(hist);
	free(x);
}

void
dump_tdist(int df)
{
	int		i, n = 100;
	float	x, val;

	for (i = 0; i < n; i++) {
		x = (float)(i - n/2) * 5 / n;
		val = Num_t2p(x, df);	// ## this is not t-dist, but t-test ###
		printf("%f %f\n", x, val);
	}
}

// BOLD response curve
RecImage *genResp(RecLoop *lp, float dt)
{
	RecImage	*img;
	int			i, n;
	float		x, a;
	float		*p;
	// model
	float		pd = 1.0;
	float		d1 = 2.5;	// sec
	float		d2 = 3.5;	// sec
	float		t1 = 5.0;	// sec
	float		t2 = 2.0;	// sec

	img = [RecImage imageOfType:RECIMAGE_REAL withLoops:lp, nil];
	p = [img data];
	n = [lp dataLength];
	for (i = 0; i < n; i++) {	// x 10msec
		x = dt * i;	// sec
		if (x >= 0 && x < d1) {
			a = 0;
		} else
		if (x < 10 + d2) {
			a = 1.0 - exp(-(x - d1) / t1);
		} else {
			a = 1.0 - exp(-(10.0 + d2 - d1) / t1);
			a *= exp(-(x - 10 - d2) / t2);
		}
		a = pd * (a - 0.4);
		p[i] = a;
	//	printf("%f %f\n", x, a);
	//	printf("%d %f\n", i, a);
	}
	// ### not done yet
	
	return img;
}

