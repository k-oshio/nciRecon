//
//	phase-gradient / k-edge
//

#import <RecKit/RecKit.h>


RecImage	*pGrad(RecImage *f, RecLoop *lp);
RecImage	*kEdge(RecImage *f, RecLoop *lp);
RecImage	*kSlit(RecImage *f, RecLoop *lp);
RecImage	*kHalf(RecImage *f, RecLoop *lp);
RecImage	*kSpaceEn(RecImage *f, RecLoop *lp, int k);
RecImage	*kSpaceRot(RecImage *f, RecLoop *lp, float rot);
RecImage	*pDif(RecImage *f, RecLoop *lp);

RecImage	*kEdge1d(RecImage *f, int xFlag);
RecImage	*kEdge1d_2(RecImage *f, int xFlag);
RecImage	*kFiltW(RecImage *f);

void		test0(RecImage *f);	// original pGrad (Liu method)
void		test1(RecImage *f);	// single, x/y/sum
void		test2(RecImage *f);	// shift spectrum
void		test3(RecImage *f);	// slope
void		test4(RecImage *f);	// take real
void		test5(RecImage *f);	// 1D k-edge
void		test6(RecImage *f);	// 1D k-edge 2
void		test7(RecImage *f);	// k0 energy ? -> X
void		test8(RecImage *f);	// phase difference
void		test9(RecImage *f);	// k-space energy
void		test10(void);		// PSF
void		test11(RecImage *f);	// histogram based
void		test12(RecImage *f);	// rot spectrum
void		test13(RecImage *f);	// dipole model

void		step_filt(RecImage *f, int xFlag);
void		state_filt(RecImage *f, int xFlag);
void		zerofill(RecImage *f, int m);
void		sctPlot(RecImage *f1, RecImage *neg);
int			sign = 1;	// f2 - f1

int
main(int ac, char *av[])
{
	RecImage	*f;
	char		fname[256];
	NSString	*path;

	if (ac > 1) {
		path = [NSString stringWithUTF8String:av[1]];
		f = [RecImage imageWithKOImage:path];
	} else {	
//		f = [RecImage imageWithKOImage:path];
//		f = [RecImage imageWithKOImage:@"raw_gr.block"];
//		f = [RecImage imageWithKOImage:@"PH_k.img"];
//		f = [RecImage imageWithKOImage:@"raw0420.block"];
		f = [RecImage imageWithKOImage:@"raw0508_1.block"];
//		f = [RecImage imageWithKOImage:@"raw0508_3.block"];
//		f = [RecImage imageWithKOImage:@"raw_imag.img"];
//		f = [RecImage imageWithKOImage:@"raw0625.block"];
	}
	if (f == nil) {
		printf("Couldn't open image file [%s]\n", fname);
		exit(0);
	}

//	test0(f);	// Liu
//	test1(f);	// single shift
//	test2(f);	// shift spectrum
//	test3(f);	// slope
//	test4(f);	// make pure real image
//	test5(f);	// 1D (various k-space energy modes)
//	test6(f);	// 1D (more modes ### move hist here)
//	test7(f);	// 2D k0
//	test8(f);	// phase difference
//	test9(f);	// simplified k-space energy
//	test10();	// PSF
//	test11(f);	// histogram based method
	test12(f);	// rot spectrum
//	test13(f);	// dipole model

	return 0;
}

// phase grad (Liu method)
void
test0(RecImage *f)
{
	RecImage	*fx, *fy;
	RecLoop		*kx, *ky;
//printf("Liu method\n");
	kx = [f xLoop];
	ky = [f yLoop];

	fx = pGrad(f, kx);
	[fx saveAsKOImage:@"IMG_x"];
	fy = pGrad(f, ky);
	[fy saveAsKOImage:@"IMG_y"];
}

// single shift
void
test1(RecImage *f)
{
	RecImage	*fx, *fy, *g;
	RecLoop		*kx, *ky;
	int			m = 2;

	kx = [f xLoop];
	ky = [f yLoop];

	g = [f copy];

	fx = kEdge(f, kx);
	[fx saveAsKOImage:@"IMG_x"];
exit(0);	// testing kEdge()
	fy = kEdge(f, ky);
	[fy saveAsKOImage:@"IMG_y"];
	[fx maxImage:fy];
	[fx saveAsKOImage:@"IMG_min"];

	[fx fft2d:REC_INVERSE];
	[fx fermi];

	// zerofill
	zerofill(fx, 2);
	[fx fft2d:REC_FORWARD];
	[fx saveAsKOImage:@"IMG_min_hires"];

	[g fermi];
	zerofill(g, 2);
	[g fft2d:REC_FORWARD];
	[g saveAsKOImage:@"IMG_mag_hires"];
}

void
zerofill(RecImage *f, int m)
{
	int			newDim;
	RecLoop		*newLoop;
	RecLoop		*kx = [f xLoop];
	RecLoop		*ky = [f yLoop];
	
	newDim = [ky dataLength] * m;
	newLoop = [RecLoop loopWithDataLength:newDim];
	[f replaceLoop:ky withLoop:newLoop];
	newDim = [kx dataLength] * m;
	newLoop = [RecLoop loopWithDataLength:newDim];
	[f replaceLoop:kx withLoop:newLoop];
}

// shift spectrum
void
test2(RecImage *f)
{
	RecImage		*fs, *fx, *fy, *ff;
	RecLoop			*kx, *ky, *sft;
	RecLoopControl	*lc;
	RecLoopIndex	*li;
	RecLoopState	*ls;
	NSString		*path;
	int				i, n = 20;
	float			d, dstep = 0.2;

	ff = [f copy];
	[ff fft2d:REC_FORWARD];
	[ff saveAsKOImage:@"IMG_mag"];

	kx = [f xLoop];
	ky = [f yLoop];

	sft = [RecLoop loopWithName:@"sft" dataLength:n+1];
	fx = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:sft, ky, kx, nil];
	fy = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:sft, ky, kx, nil];
	lc = [fx control];
	li = [lc loopIndexAtIndex:0];
	[li setActive:NO];
	ls = [li state];
	[lc rewind];
	for (i = 0; i <= n; i++) {
		[ls setCurrent:i];
		d = (i - n/2) * dstep + 0.5;
		fs = [f copy];
		[fs rotate1d:kx by:d];
		[fs rotate1d:ky by:d];
		ff = kEdge(fs, kx);
		[fx copyImage:ff withControl:lc];
		fs = [f copy];
		[fs rotate1d:kx by:d];
		[fs rotate1d:ky by:d];
		ff = kEdge(fs, ky);
		[fy copyImage:ff withControl:lc];
	}
	[fx saveAsKOImage:@"IMG_x_shift"];
	[fy saveAsKOImage:@"IMG_y_shift"];
	[fx minImage:fy];
	[fx saveAsKOImage:@"IMG_min_shift"];
}

// slope
void
test3(RecImage *f)
{
	RecImage	*fx, *fy, *ff;
	RecImage	*fp, *fm;
	RecImage	*mask, *maskh;
	RecLoop		*kx, *ky;
	float		d = 0.30;
	float		thres = 0.01;

	kx = [f xLoop];
	ky = [f yLoop];
	mask = [f copy];
	[mask fft2d:REC_FORWARD];
	[mask saveAsKOImage:@"IMG_mag"];
	[mask thresAt:thres];

	fp = [f copy];
	[fp rotate1d:kx by:0.5 - d];
	[fp rotate1d:ky by:0.5 - d];
	fm = [f copy];
	[fm rotate1d:kx by:0.5 + d];
	[fm rotate1d:ky by:0.5 + d];

	ff = kEdge(fm, kx);
	fx = kEdge(fp, kx);
	[fx subImage:ff];
	[fx saveAsKOImage:@"IMG_x"];
	ff = kEdge(fm, ky);
	fy = kEdge(fp, ky);
	[fy subImage:ff];
	[fy saveAsKOImage:@"IMG_y"];
	[fx maxImage:fy];

//	[fx scaleByImage:mask];
	[fx saveAsKOImage:@"IMG_max"];

	[fx fft2d:REC_INVERSE];
	[fx fermi];
//[fx saveAsKOImage:@"IMG_min_fermi"];

	// mag hires
	mask = [f copy];
	[mask fermi];
	zerofill(mask, 2);
	[mask fft2d:REC_FORWARD];
	[mask saveAsKOImage:@"IMG_mag_hires"];
	[mask thresAt:thres];

	// min hires
	zerofill(fx, 2);
	[fx fft2d:REC_FORWARD];
	[fx scaleByImage:mask];
	[fx saveAsKOImage:@"IMG_max_hires"];
}

void
test4(RecImage *f)	// make real image for testing
{
	[f fft2d:REC_FORWARD];
	[f takeRealPart];
//	[f magnitude];

//	[f makeComplex];
	[f fft2d:REC_INVERSE];
//	[f saveAsKOImage:@"raw_real.img"];
	[f saveAsKOImage:@"raw_imag.img"];
}

// 1D version (various modes)
void
test5(RecImage *f)
{
	RecImage	*fx, *fy, *fm;

//	[f rotate1d:[f xLoop] by:-0.5];
//	[f rotate1d:[f yLoop] by:-0.5];

	fm = [f copy];
	[fm fft2d:REC_FORWARD];
	[fm saveAsKOImage:@"IMG_mag"];

	fx = kEdge1d(f, YES);
	[fx saveAsKOImage:@"IMG_x"];
	fy = kEdge1d(f, NO);
	[fy saveAsKOImage:@"IMG_y"];

	[fx maxImage:fy];
	[fx saveAsKOImage:@"IMG_max"];
}

// 1D version (various modes)
void
test6(RecImage *f)
{
	RecImage	*fx, *fy, *fm;

	fm = [f copy];
	[fm fft2d:REC_FORWARD];
	[fm saveAsKOImage:@"IMG_mag"];

	fx = [f copy];
	step_filt(fx, YES);
	[fx saveAsKOImage:@"IMG_f1x"];
//	[fx subImage:fm];
	[fm addConst:1000.0];
	[fx divImage:fm];
	[fx saveAsKOImage:@"IMG_f1x_sub"];

// 2D filt
	state_filt(fx, YES);
	[fx saveAsKOImage:@"IMG_x"];

	fy = [f copy];
	step_filt(fy, NO);
	[fy saveAsKOImage:@"IMG_f1y"];
// 2D filt
	state_filt(fy, NO);
	[fy saveAsKOImage:@"IMG_y"];

// scatter plot

}

// 2D k-0
void
test7(RecImage *f)
{
	RecImage		*f1, *f2;
	RecImage		*pos, *neg;
	RecLoop			*kx, *ky;
	int				i, k;
	int				len, pixSize, dataLength;
	int				xDim, yDim;
	float			*p, *q;
	RecLoopControl	*lc;

	f1 = [f copy];
	f2 = [f copy];
	lc = [f1 control];
	[lc deactivateXY];

	pixSize = [f1 pixSize];
	dataLength = [f1 dataLength];
	xDim = [f xDim];
	yDim = [f yDim];

	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f1 currentDataWithControl:lc];
		q = p + dataLength;
		p[(yDim / 2) * xDim + xDim / 2] = 0;
		q[(yDim / 2) * xDim + xDim / 2] = 0;
		[lc increment];
	}

	[f1 saveAsKOImage:@"IMG_f1"];
	[f2 saveAsKOImage:@"IMG_f2"];

	[f1 fft2d:REC_FORWARD];
	[f2 fft2d:REC_FORWARD];
	[f1 magnitude];
	[f2 magnitude];
[f1 saveAsKOImage:@"IMG_f1_mg"];

// === f1 + f2 ===
	[f1 subImage:f2];
[f1 saveAsKOImage:@"IMG_sub_mg"];
	[f1 divImage:f2];
	[f1 makeComplex];

	[f1 saveAsKOImage:@"IMG_xy"];
}

void
test8(RecImage *f)
{
	RecImage	*fx, *fy;
	RecLoop		*kx, *ky;

	[f fft2d:REC_FORWARD];

	kx = [f xLoop];
	ky = [f yLoop];

	fx = pDif(f, kx);
	[fx saveAsKOImage:@"IMG_x"];
	fy = pDif(f, ky);
	[fy saveAsKOImage:@"IMG_y"];
}

// simplified k-space energy
void
test9(RecImage *f)
{
	RecImage		*fs, *fx, *fy, *ff;
	RecLoop			*kx, *ky, *sft;
	RecLoopControl	*lc;
	RecLoopState	*sft_st, *sl_st;
	NSString		*path;
	int				i, len, n = 100;
	float			r;
	int				index = 1;	// image#

	ff = [f copy];
	[ff fft2d:REC_FORWARD];
	[ff saveAsKOImage:@"IMG_mag"];

	kx = [f xLoop];
	ky = [f yLoop];
	len = [kx dataLength];
	if (n > len) {
		n = len;
	}

	sft = [RecLoop loopWithName:@"sft" dataLength:n];
	fx = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:sft, ky, kx, nil];
	fy = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:sft, ky, kx, nil];
	fs = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ky, kx, nil];
	lc = [f control];	// slice, ky, kx
	[lc insertLoop:sft atIndex:0];				// sft, slice, ky, kx
	sft_st = [lc loopStateAtIndex:0];
	sl_st = [lc loopStateAtIndex:1];
	[[lc loopIndexAtIndex:0] setActive:NO];
	[[lc loopIndexAtIndex:1] setActive:NO];

	[sl_st setCurrent:index];		// select slice
	[lc rewind];
	for (i = 0; i < n; i++) {
		[sft_st setCurrent:i];	// set shift
		r = (float)(i - n/2) * 0.1;
	//	printf("%d\r", i); fflush(stdout);
		// x
		[fs copyImage:f withControl:lc];
		ff = kSpaceRot(fs, kx, r);
		[fx copyImage:ff withControl:lc];
		// y
		[fs copyImage:f withControl:lc];
		ff = kSpaceRot(fs, ky, r);
		[fy copyImage:ff withControl:lc];
	}
	[fx saveAsKOImage:@"IMG_x_kse"];
	[fy saveAsKOImage:@"IMG_y_kse"];
}

void
test10()
{
	RecImage	*f;
	RecLoop		*ky, *kx;
	float		*p, *q, x;
	int			xdim, ydim, skip;
	int			i, j;
	int			mode = 2;	// 0:const, 1:sinc 2:pgrad

	xdim = ydim = 256;
	kx = [RecLoop loopWithDataLength:xdim];
	ky = [RecLoop loopWithDataLength:ydim];
	f = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:ky, kx, nil];

	p = [f data];
	q = p + [f dataLength];
	skip = [f skipSizeForLoop:ky];

	for (i = 0; i < ydim; i++) {
		for (j = 0; j < xdim; j++) {
			if (j > i) {
				switch (mode) {
				// const
				case 0 :
					p[j] = 1.0;
					break;
				// sinc
				case 1 :
					x = M_PI * 2 * (j - xdim/2) / xdim;
					if (x != 0) {
						p[j] = sin(x) / x;
					} else {
						p[j] = 1.0;
					}
					break;
				case 2 :
					x = M_PI * 2 * (j - xdim/2) / xdim + 1.0;
					if (x != 0) {
						p[j] = sin(x) / x;
					} else {
						p[j] = 1.0;
					}
					break;
				}
			} else {
				p[j] = 0.0;
			}
			q[j] = 0.0;
		}
		p += skip;
		q += skip;
	}
	[f saveAsKOImage:@"IMG_PSF1"];

	[f fft1d:kx direction:REC_FORWARD];
	[f saveAsKOImage:@"IMG_PSF2"];
}

// ####
void
test11(RecImage *f)	// histogram based
{
	RecImage	*fx, *fy, *fm;
	RecLoop		*kx, *ky;

	fm = [f copy];
	[fm fft2d:REC_FORWARD];	// mag


	fx = pDif(f, kx);
	[fx saveAsKOImage:@"IMG_x"];
	fy = pDif(f, ky);
	[fy saveAsKOImage:@"IMG_y"];
}

// rot spectrum
// ### extract a slice (stim phase)
// ### subtract or div const part (mag ?)
void
test12(RecImage *f)
{
	RecImage		*fs, *fx, *ff;
	RecImage		*rotMap, *rotParam;
	RecLoop			*kx, *ky, *rot;
	int				i, n = 360;
	float			th, *p;

	ff = [f copy];
	[ff fft2d:REC_FORWARD];
	[ff saveAsKOImage:@"IMG_mag"];

	kx = [f xLoop];
	ky = [f yLoop];

	rot = [RecLoop loopWithName:@"rot" dataLength:n];
	rotParam = [RecImage imageOfType:RECIMAGE_REAL withLoops:rot, nil];

	fx = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:rot, ky, kx, nil];
	p = [rotParam data];
	for (i = 0; i < n; i++) {
		p[i] = (float)i * M_PI * 2 / n;
	}
	rotMap = [fx mapForRotate:rotParam];
	[fx resample:f withMap:rotMap];

	ff = kEdge(fx, kx);
//	ff = kSlit(fx, kx);
//	ff = kHalf(fx, kx);
//    ff = kFiltW(fx);
	[ff saveAsKOImage:@"IMG_flt"];
	p = [rotParam data];
	for (i = 0; i < n; i++) {
		p[i] = -(float)i * M_PI * 2 / n;
	}
	rotMap = [fx mapForRotate:rotParam];
	[fx resample:ff withMap:rotMap];
	[fx saveAsKOImage:@"IMG_rot"];

//	fs = [fx maxForLoop:rot];
//	[fs saveAsKOImage:@"IMG_rot_max"];
//	[fx magnitudeSq];
//	[fx magnitudeSq];
//	fs = [fx sumForLoop:rot];
//	[fs saveAsKOImage:@"IMG_rot_sum"];
}

// dipole model
void
test13(RecImage *f)
{
	RecImage		*f1, *f2;
	RecImage		*mask, *mag;
	RecLoopControl	*lc;
	int				xDim = [f xDim];
	int				yDim = [f yDim];
	int				i, j, k, ix;
	int				n, dataLength;
	float			*p, *q, *pp, mg, mx;
	float			x, y, r, cs, w, invr, invi, re, im;
	int				phsmode = 0;	// 0: filter imag, then take phase 1: take phs, then filter
	int				mode = 2;		// 0: cos, 1: 1/cos 2: 1/(cos + s) 3: r^2 cos / (cos^2 + s^2 r^2)
	float			mgs = 0.1;
	float			sig = 0.2;		// noise in imag part [0..1]
	float			sd = 1.0;		// PSF width
	int				xFlag = NO;

	dataLength = [f dataLength];
	lc = [f control];
	[lc deactivateXY];
	n = [lc loopLength];

// mask image
	mag = [f copy];
	[mag fft2d:REC_FORWARD];
	[mag magnitude];
	mx = [mag maxVal];
	mgs *= mx;
	[mag saveAsKOImage:@"IMG_mag"];
	mask = [mag copy];
	[mask thresAt:0.05];
	[mask saveAsKOImage:@"IMG_mask"];

/* equivalent to taking imag part
	f1 = [f copy];
	[f1 rotate:2];
	[f1 rotate1d:[f xLoop] by:-1.0];	// center is n/2
	[f1 rotate1d:[f yLoop] by:-1.0];	// center is n/2
	[f1 conj];
	[f1 negate];
	[f1 addImage:f];
*/

/* take imag part */
	if (phsmode == 0) {
		f1 = [f copy];
		[f1 fft2d:REC_FORWARD];
		p = [f1 data];
		for (i = 0; i < dataLength; i++) {
			p[i] = 0;
		}
		[f1 fft2d:REC_INVERSE];
		[f1 saveAsKOImage:@"IMG_imag_ft"];
	} else
/* take phase, and mask */
	if (phsmode == 1) {
		f1 = [f copy];
		[f1 fft2d:REC_FORWARD];
		p = [f1 data];
		q = p + dataLength;
		pp = [mag data];
		for (i = 0; i < dataLength; i++) {
			mg = pp[i] + mgs;
			p[i] = 0;
			q[i] /= mg;
		}
		[f1 saveAsKOImage:@"IMG_imag_tmp"];
		[f1 fft2d:REC_INVERSE];
		[f1 saveAsKOImage:@"IMG_imag_ft"];
	}

// PSF filter
	[lc rewind];
	for (k = 0; k < n; k++) {
		p = [f1 currentDataWithControl:lc];
		q = p + dataLength;
		for (i = ix = 0; i < yDim; i++) {
			y = i - yDim/2;
			for (j = 0; j < xDim; j++, ix++) {
				x = j - xDim/2;
				r = sqrt(x*x + y*y);
				if (r < 0.2) {
					w = 0;
				} else {
				// inverse filter
					if (xFlag) {
						cs = x/r;
					} else {
						cs = y/r;
					}
					w = exp(-r*r/(sd*xDim*yDim));	// -> sinc
					switch (mode) {
					case 0:
						invr = cs;
						p[ix] *= w * invr;
						q[ix] *= w * invr;
						break;
					case 1:
						if (cs > 0.1) {
							invr = 1/cs;
						} else {
							invr = 0;
						}
						p[ix] *= w * invr;
						q[ix] *= w * invr;
						break;
					case 2:
						invr = cs / (cs*cs + sig*sig);
						p[ix] *= w * invr;
						q[ix] *= w * invr;
						break;
					case 3:
						invr = r*r * cs / (cs*cs + sig*sig*r*r);
						p[ix] *= w * invr;
						q[ix] *= w * invr;
						break;
					case 4:
						invr = cs * w / (cs*cs*w*w + sig*sig);
						p[ix] *= invr;
						q[ix] *= invr;
						break;
					}
				}
			}
		}
		[lc increment];
	}
	[f1 saveAsKOImage:@"IMG_imag_flt"];	// change name ###

// FT
	[f1 fft2d:REC_FORWARD];
	if (phsmode == 1) {
	//	[f1 scaleByImage:mask];
	}
	[f1 saveAsKOImage:@"IMG_imag"];	// change name ###

/* take phase */
	if (phsmode == 0) {
		[lc rewind];
		for (k = 0; k < n; k++) {
			pp = [mag currentDataWithControl:lc];
			p = [f1 currentDataWithControl:lc];
			q = p + dataLength;
			for (i = 0; i < xDim * yDim; i++) {
				mg = pp[i] + mgs;
				p[i] /= mg;
				q[i] /= mg;
			}
			[lc increment];
		}
	//	[f1 limitToVal:100000.0];
		[f1 saveAsKOImage:@"IMG_tmp_ft2"];	// change name ###
	}

	f1 = [f1 oversample];
	[f1 saveAsKOImage:@"IMG_tmp_ft2_hires"];	// change name ###
	mag = [mag oversample];
	[mag saveAsKOImage:@"IMG_mag_hires"];	// change name ###
}

RecImage *
pGrad(RecImage *f, RecLoop *lp)
{
	RecImage		*f1, *f2;
	RecImage		*pos, *neg;
	RecLoop			*kx, *ky;
	int				i, k;
	int				len, pixSize, dataLength;
	float			sft = 0.4;
	int				skip;
	float			*p;
	RecLoopControl	*lc;

//	[f rotate1d:[f xLoop] by:0.5];
//	[f rotate1d:[f yLoop] by:0.5];
	f1 = [f copy];
	f2 = [f copy];
	lc = [RecLoopControl outerLoopControlForImage:f1 forLoop:lp];

	len = [lp dataLength];
	skip = [f1 skipSizeForLoop:lp];
	pixSize = [f1 pixSize];
	dataLength = [f1 dataLength];

	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f1 currentDataWithControl:lc];
		for (k = 0; k < pixSize; k++) {
			vDSP_vclr(p + skip*(int)(len*(1.0 - sft)), skip, (int)(len * sft));
			p += dataLength;
		}
		[lc increment];
	}
	[f1 saveAsKOImage:@"IMG_f1"];

	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f2 currentDataWithControl:lc];
		for (k = 0; k < pixSize; k++) {
			vDSP_vclr(p, skip,  (int)(len * sft));
			p += dataLength;
		}
		[lc increment];
	}
	[f2 saveAsKOImage:@"IMG_f2"];

kx = [RecLoop loopWithDataLength:256];
ky = [RecLoop loopWithDataLength:256];
[kx setDataLength:256];
[ky setDataLength:256];
[f1 replaceLoop:[f1 xLoop] withLoop:kx];
[f1 replaceLoop:[f1 yLoop] withLoop:ky];
[f2 replaceLoop:[f2 xLoop] withLoop:kx];
[f2 replaceLoop:[f2 yLoop] withLoop:ky];
	[f1 fft2d:REC_FORWARD];
	[f2 fft2d:REC_FORWARD];
	[f1 magnitude];
	[f2 magnitude];

// === f1 + f2 ===
	pos = [f1 copy];
	[pos addImage:f2];

	neg = [f2 copy];
	[neg subImage:f1];


//[pos saveAsKOImage:@"IMG_pos"];
[neg saveAsKOImage:@"IMG_neg"];

	[neg divImage:pos];
	[neg makeComplex];

	return neg;
}

// kEdge
RecImage *
kEdge(RecImage *f, RecLoop *lp)
{
	RecImage		*f1, *f2;
	RecImage		*pos, *neg;
	int				i, k;
	int				len, pixSize, dataLength;
	int				skip;
	float			*p;
	RecLoopControl	*lc;
	f1 = [f copy];
	f2 = [f copy];
	lc = [RecLoopControl outerLoopControlForImage:f1 forLoop:lp];

	len = [lp dataLength];
	skip = [f1 skipSizeForLoop:lp];
	pixSize = [f1 pixSize];
	dataLength = [f1 dataLength];

	// f1
	f1 = [f copy];
	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f1 currentDataWithControl:lc];
		for (k = 0; k < pixSize; k++) {
			// L
			//vDSP_vclr(p + skip*(len/2), skip, len/2);
			// R-1
			vDSP_vclr(p, skip, len/2 + 1);
			p += dataLength;
		}
		[lc increment];
	}

	// f2
	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f2 currentDataWithControl:lc];
		for (k = 0; k < pixSize; k++) {
			// R
			vDSP_vclr(p, skip, len/2);
			// L+1
			//vDSP_vclr(p + skip*(len/2) + 1, skip, len/2 - 1);
			p += dataLength;
		}
		[lc increment];
	}
//	[f1 saveAsKOImage:@"IMG_f1"];
//	[f2 saveAsKOImage:@"IMG_f2"];

	[f1 fft2d:REC_FORWARD];
	[f2 fft2d:REC_FORWARD];
	[f1 saveAsKOImage:@"IMG_f1_ft"];
	[f2 saveAsKOImage:@"IMG_f2_ft"];
	[f1 magnitude];		// default
	[f2 magnitude];		// default

// === f1 only version
    [f1 addConst:10000.0];
    [f1 invert];
    [f1 makeComplex];
[f1 saveAsKOImage:@"IMG_f1_inv"];
    return f1;

// === f1 + f2 (probably this is correct) ===
	pos = [f1 copy];
	[pos addImage:f2];

// === mag (not this one) ====
//	pos = [f copy];
//	[pos fft2d:REC_FORWARD];
//	[pos magnitude];

//	neg = [f2 copy];
//	[neg subImage:f1];
	neg = [f2 copy];
	
	[pos saveAsKOImage:@"IMG_pos"];
//	[neg saveAsKOImage:@"IMG_neg"];

	[neg divImage:pos];
	[neg addConst:-0.5];

	[neg makeComplex];

	return neg;
}

// kSlit
RecImage *
kSlit(RecImage *f, RecLoop *lp)
{
	RecImage		*f1, *f2;
	RecImage		*pos, *neg;
	int				i, k;
	int				len, pixSize, dataLength;
	int				skip;
	float			*p;
	RecLoopControl	*lc;
	f1 = [f copy];
	f2 = [f copy];
	lc = [RecLoopControl outerLoopControlForImage:f1 forLoop:lp];

	len = [lp dataLength];
	skip = [f1 skipSizeForLoop:lp];
	pixSize = [f1 pixSize];
	dataLength = [f1 dataLength];

	// f1 (slit)
	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f1 currentDataWithControl:lc];
		for (k = 0; k < pixSize; k++) {
			p[(len/2)*skip] = 0;
			p += dataLength;
		}
		[lc increment];
	}

	// f2 (mag), do nothing
	[f1 saveAsKOImage:@"IMG_f1"];
	[f2 saveAsKOImage:@"IMG_f2"];

	[f1 fft2d:REC_FORWARD];
	[f2 fft2d:REC_FORWARD];
	[f1 saveAsKOImage:@"IMG_f1_ft"];
	[f2 saveAsKOImage:@"IMG_f2_ft"];
	[f1 magnitude];		// default
	[f2 magnitude];		// default

// === f1 + f2 (probably this is correct) ===
	pos = [f1 copy];
	[pos addImage:f2];

// === mag (not this one) ====
//	pos = [f copy];
//	[pos fft2d:REC_FORWARD];
//	[pos magnitude];

	neg = [f2 copy];
	[neg subImage:f1];
	
	[pos saveAsKOImage:@"IMG_pos"];
	[neg saveAsKOImage:@"IMG_neg"];

	[neg divImage:pos];
	[neg makeComplex];

	return neg;
}

// kSlit
RecImage *
kHalf(RecImage *f, RecLoop *lp)
{
	RecImage		*f1, *f2;
	RecImage		*pos, *neg;
	int				i, k;
	int				len, pixSize, dataLength;
	int				skip;
	float			*p;
	RecLoopControl	*lc;
	f1 = [f copy];
	f2 = [f copy];
	lc = [RecLoopControl outerLoopControlForImage:f1 forLoop:lp];

	len = [lp dataLength];
	skip = [f1 skipSizeForLoop:lp];
	pixSize = [f1 pixSize];
	dataLength = [f1 dataLength];

	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f1 currentDataWithControl:lc];
		for (k = 0; k < pixSize; k++) {
			vDSP_vclr(p, skip, len/2);
			p[(len/2)*skip] *= 0.5;
			p += dataLength;
		}
		[lc increment];
	}
	[f1 fft2d:REC_FORWARD];
	[f1 saveAsKOImage:@"IMG_f0_ft"];

	[f1 takeRealPart];
	[f2 fft2d:REC_FORWARD];
	[f2 takeRealPart];
//	[f2 multByConst:0.5];
//	[f1 subImage:f2];
	[f2 magnitude];
	[f2 addConst:100.0];
	[f1 divImage:f2];
	[f1 makeComplex];

	return f1;
}

// kEdge 1D version
RecImage *
kEdge1d(RecImage *f, int xFlag)
{
	RecImage		*f1, *f2, *ff;
	RecImage		*pos, *neg;
	int				i, j, k, ix;
	int				len, pixSize, dataLength;
	int				skip, w;
	float			*p, *q;
	RecLoopControl	*lc;
	RecLoop			*lp, *other;
	int				mode = 1;
		//	0: k-edge (2 side), 1: 1 side (default), 2: k-0, 3: lo-pass
		//	4:(removed) 5: (removed), 6: Liu's pGrad 7: k-edge + level correction
		//	8: f1 only
	float			lv = 0.5;
	float			r1 = 1.0;
	float			r2 = 1.0;

	f1 = [f copy];
	f2 = [f copy];

	if (xFlag) {
		lp = [f xLoop];
		other = [f yLoop];
	} else {
		lp = [f yLoop];
		other = [f xLoop];
	}
	lc = [RecLoopControl outerLoopControlForImage:f1 forLoop:lp];

	len = [lp dataLength];
	skip = [f1 skipSizeForLoop:lp];
	pixSize = [f1 pixSize];
	dataLength = [f1 dataLength];

// f1
	[f1 fft1d:other direction:REC_FORWARD];
	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f1 currentDataWithControl:lc];
		q = p + dataLength;
		switch (mode) {
		case 0 :	// k-edge
//		case 5 :	// imag-energy
			vDSP_vclr(p + skip*(len/2), skip, len/2);
			vDSP_vclr(q + skip*(len/2), skip, len/2);
			break;
		case 7 :
			lv = -0.4;
			vDSP_vclr(p, skip, len/2);
			vDSP_vclr(q, skip, len/2);
			p[len/2 * skip] *= lv;
			q[len/2 * skip] *= lv;
			break;
		case 1 : // k-space energy dif (probably highest SNR, theoretically simpler)
		case 8 :
			vDSP_vclr(p, skip, len/2 + 1);
			vDSP_vclr(q, skip, len/2 + 1);
			break;
		case 2 : // k-0
			p[len/2 * skip] *= lv;
			q[len/2 * skip] *= lv;
			break;
		case 3 : // test high freq
			w = 20;
			vDSP_vclr(p, skip, w);
			vDSP_vclr(q, skip, w);
			break;
		case 6 :
			w = len/2 - 16;
			vDSP_vclr(p + skip*(len - w), skip, w);
			vDSP_vclr(q + skip*(len - w), skip, w);
			break;
		}


		[lc increment];
	}

// f2
	[f2 fft1d:other direction:REC_FORWARD];
	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f2 currentDataWithControl:lc];
		q = p + dataLength;
		switch (mode) {
		case 0 :	// k-edge
		case 1 :	// k-space energy dif
//		case 5 :	// imag/power
		case 7 :
		case 8 :
			vDSP_vclr(p, skip, len/2);
			vDSP_vclr(q, skip, len/2);
			break;
		case 2 :
			p[(len/2 - 1) * skip] *= lv;
			q[(len/2 - 1) * skip] *= lv;
			break;
		case 3 : // test high freq
			w = 20;
			vDSP_vclr(p + (len - w) * skip, skip, w);
			vDSP_vclr(q + (len - w) * skip, skip, w);
			break;
		case 6 :
			w = len/2 - 16;
			vDSP_vclr(p, skip, w);
			vDSP_vclr(q, skip, w);
			break;
		}
		[lc increment];
	}
	[f1 saveAsKOImage:@"IMG_f1"];
	[f2 saveAsKOImage:@"IMG_f2"];

	[f1 fft1d:lp direction:REC_FORWARD];
	[f2 fft1d:lp direction:REC_FORWARD];

	[f1 saveAsKOImage:@"IMG_f1_ft"];
	[f2 saveAsKOImage:@"IMG_f2_ft"];
	ff = [f2 copy];
	switch (mode) {
	default :
		[f1 magnitude];
		[f2 magnitude];
		break;
	}
	[f1 saveAsKOImage:@"IMG_f1_mg"];
	[f2 saveAsKOImage:@"IMG_f2_mg"];


//======================== refine this part to compare combination of f1 & f2 =====
	// === pos ===
	switch (mode) {
	case -1 :
	// === mag ====
		pos = [f copy];
		[pos fft2d:REC_FORWARD];
		[pos magnitude];
		break;
	case 9 :
	// === exp ====
		pos = [f2 copy];
		[pos multByConst:1/r1];
	//	[pos exp];
		break;
	default :
	// === f1 + f2 ====
		pos = [f1 copy];
		[pos addImage:f2];
		break;
	}
	// === neg ===
	switch (mode) {
	case 8 :
		neg = [f2 copy];
		[neg divImage:pos];
		break;
	case 9 :
	// === exp ====
		neg = [f1 copy];
		[neg multByConst:-1/r2];
	//	[neg exp];
		[neg negate];
		[neg addConst:1.0];
		break;
	default :
		neg = [f2 copy];
		[neg subImage:f1];
		[neg divImage:pos];
		break;
	}
	[pos saveAsKOImage:@"IMG_pos"];
	[neg saveAsKOImage:@"IMG_neg"];
//======================== refine this part to compare combination of f1 & f2 =====

	[neg makeComplex];

// 2D scatter plot
	if (xFlag) sctPlot(ff, neg);
	return neg;
}

// inter-pixel phase difference
RecImage *
pDif(RecImage *f, RecLoop *lp)
{
	RecImage		*df;
	int				i, j, k;
	int				len, pixSize, dataLength;
	int				skip;
	float			*p, *q;
	float			*pp, *qq, re1, re2, im1, im2;
	RecLoopControl	*lc;

	df = [f copy];
	lc = [RecLoopControl outerLoopControlForImage:f forLoop:lp];

	len = [lp dataLength];
	skip = [f skipSizeForLoop:lp];
printf("skip = %d\n", skip);
	pixSize = [f pixSize];
	dataLength = [f dataLength];

	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f currentDataWithControl:lc];
		q = p + dataLength;
		pp = [df currentDataWithControl:lc];
		qq = pp + dataLength;
		pp[0] = qq[0] = 0;
		for (j = 1; j < len; j++) {
			k = j * skip;
			re1 = p[k - skip];
			re2 = p[k];
			im1 = q[k - skip];
			im2 = q[k];
			pp[k] = re1 * re2 + im1 * im2;
			qq[k] = im1 * re2 - re1 * im2;
			pp[k] = atan2(qq[k], pp[k]);
			qq[k] = 0;
		}
		[lc increment];
	}

	return df;
}

RecImage *
kSpaceEn(RecImage *f, RecLoop *lp, int kpos)
{
	RecImage		*f1, *mg;
	int				i, k;
	int				len, pixSize, dataLength;
	int				skip;
	float			*p;
	RecLoopControl	*lc;

	f1 = [f copy];
//	mg = [f copy];
//	[mg magnitude];
	lc = [RecLoopControl outerLoopControlForImage:f1 forLoop:lp];

	len = [lp dataLength];
	skip = [f1 skipSizeForLoop:lp];
	pixSize = [f1 pixSize];
	dataLength = [f1 dataLength];

	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f1 currentDataWithControl:lc];
		for (k = 0; k < pixSize; k++) {
			vDSP_vclr(p + skip*kpos, skip, len - kpos);
			p += dataLength;
		}
		[lc increment];
	}
//	[f1 saveAsKOImage:@"IMG_f1"];
	[f1 fft2d:REC_FORWARD];
	[f1 magnitude];

//	[f1 divImage:mg];
	[f1 makeComplex];

	return f1;
}

RecImage *
kSpaceRot(RecImage *f, RecLoop *lp, float rot)
{
	RecImage		*f1, *mg;
	int				i, j, k, ix;
	int				len, pixSize, dataLength;
	int				skip;
	float			*p;
	RecLoopControl	*lc;

	f1 = [f copy];
	lc = [RecLoopControl outerLoopControlForImage:f1 forLoop:lp];

	len = [lp dataLength];
	skip = [f1 skipSizeForLoop:lp];
	pixSize = [f1 pixSize];
	dataLength = [f1 dataLength];

	[lc rewind];
	[f1 rotate1d:lp by:-rot];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f1 currentDataWithControl:lc];
		for (k = 0; k < pixSize; k++) {
			// zero left
			//vDSP_vclr(p, skip, len/2);
			// zero right
			//vDSP_vclr(p + len/2 * skip, skip, len/2);
			// lev corr
			for (j = ix = 0; j < len/2; j++, ix += skip) {
				p[ix] *= -0.2;
			}
			p += dataLength;
		}
		[lc increment];
	}
	[f1 fft2d:REC_FORWARD];
//	[f1 magnitude];
//	[f1 makeComplex];

	return f1;
}

void
step_filt(RecImage *f, int xFlag)
{
	RecLoop			*lp, *other;
	RecLoopControl	*lc;
	int				i, len, skip, dataLength;
	float			*p, *q;
	BOOL			half_flg = YES;

	if (xFlag) {
		lp = [f xLoop];
		other = [f yLoop];
	} else {
		lp = [f yLoop];
		other = [f xLoop];
	}
	lc = [RecLoopControl outerLoopControlForImage:f forLoop:lp];

	len = [lp dataLength];
	skip = [f skipSizeForLoop:lp];
	dataLength = [f dataLength];

// f1
//	[f rotate1d:lp by:0.0];
	[f fft1d:other direction:REC_FORWARD];
	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f currentDataWithControl:lc];
		q = p + dataLength;
		vDSP_vclr(p, skip, len/2);
		vDSP_vclr(q, skip, len/2);
		if (half_flg) {
			p[len/2 * skip] *= 0.5;
			q[len/2 * skip] *= 0.5;
		}
		[lc increment];
	}
	[f fft1d:lp direction:REC_FORWARD];

	if (half_flg) {
		[f takeRealPart];
		[f multByConst:2.0];
		[f makeComplex];
	}
}

void
state_filt(RecImage *f, int xFlag)
{
	RecLoop			*lp, *other;
	RecLoopControl	*lc;
	int				i, j, ix, len, skip, dataLength;
	float			*p, *q;
	float			mp, mq;
	float			x, y, r0, r1;
	float			d0 = 0.5;
	float			d1 = 0.2;
	int				mode = 1;
	float			lv = 1.0; //1.3;

	if (xFlag) {
		lp = [f xLoop];
		other = [f yLoop];
	} else {
		lp = [f yLoop];
		other = [f xLoop];
	}
	lc = [RecLoopControl outerLoopControlForImage:f forLoop:lp];

	len = [lp dataLength];
	skip = [f skipSizeForLoop:lp];
	dataLength = [f dataLength];

// f1
	[lc rewind];
	for (i = 0; i < [lc loopLength]; i++) {
		p = [f currentDataWithControl:lc];
		q = p + dataLength;
		vDSP_sve(p, skip, &mp, len);
		vDSP_sve(q, skip, &mq, len);
		mp /= len;
		mq /= len;

		mp *= lv;
		mq *= lv;

		for (j = ix = 0; j < len; j++, ix += skip) {
			x = p[ix];
			y = q[ix];
			r0 = sqrt(x*x + y*y);
			r1 = sqrt((x - mp)*(x - mp) + (y - mq)*(y - mq));
			switch (mode) {
			case 0 :
				p[ix] = (r0 - r1) / (r0 + r1);
				q[ix] = 0;
				break;
			case 1 :
				p[ix] = r0 / (r0 + r1);
				q[ix] = 0;
				break;
			case 2 :
				p[ix] = r0*r0 / (r0*r0 + r1*r1);
				q[ix] = 0;
				break;
			case 3 :
				p[ix] = (1 + exp(-r1 / (mp * d0))) * (1 - exp(-r0 / (mp * d1)));
				q[ix] = 0;
				break;
			}
		}
		[lc increment];
	}
}

// x loop only
void
sctPlot(RecImage *f1, RecImage *neg)
{
	RecImage		*sct1, *sct2;
	RecLoop			*reLp, *imLp;
	RecLoopControl	*f1Lc, *negLc;
	RecLoopState	*f1YSt, *f1SlSt;
	RecLoopState	*negYSt, *negSlSt;
	float			*p, *q, *s, mn, mx, sm;
	float			*pp1, *qq1;
	float			*pp2, *qq2;
	int				i, j, n, len, xdim, ydim;
	int				x, y;
	int				f1Sl = 1;
	int				negSl = 1;
	int				ii, jj;

	n = 512;
	f1Lc = [f1 control];
	negLc = [neg control];
	f1SlSt = [f1Lc loopStateAtIndex:0];
	negSlSt = [negLc loopStateAtIndex:0];
	f1YSt = [f1Lc loopStateAtIndex:1];
	negYSt = [negLc loopStateAtIndex:1];
	len = [neg dataLength];
	xdim = [neg xDim];
	ydim = [neg yDim];

	reLp = 	[RecLoop loopWithDataLength:n];
	imLp = 	[RecLoop loopWithDataLength:n];
	sct1 = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:reLp, imLp, nil];
	sct2 = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:reLp, imLp, nil];
	pp1 = [sct1 data];
	qq1 = pp1 + [sct1 dataLength];
	pp2 = [sct2 data];
	qq2 = pp2 + [sct1 dataLength];

// global max
	mx = sm = 0;
	[f1SlSt setCurrent:f1Sl];
	[negSlSt setCurrent:negSl];
	p = [f1 currentDataWithControl:f1Lc];
	s = [neg currentDataWithControl:negLc];
	for (i = 0; i < xdim*ydim; i++) {
		if (mx < p[i]) mx = p[i];
		sm += s[i];
	}
	sm /= (xdim*ydim);
	[sct1 setConst:sm];
	[sct2 setConst:sm];
	for (i = 0; i < ydim; i++) {
		[f1YSt setCurrent:i];
		[negYSt setCurrent:i];
		p = [f1 currentDataWithControl:f1Lc];
		q = p + len;
		s = [neg currentDataWithControl:negLc];
		// line avg
		mn = 0;
		for (j = 0; j < xdim; j++) {
			mn += p[j];
		}
		mn /= xdim;

		for (j = 0; j < xdim; j++) {

	jj = 75;
	ii = 87;
	if ((i < ii) || (i > ii + 4) || (j < jj) || (j > jj + 4)) continue;

		// raw
			x = p[j] / mx * n + n/2;
			if (x < 0) x = 0;
			if (x >= n) x = n-1;
			y = q[j] / mx * n + n/2;
			if (y < 0) y = 0;
			if (y >= n) y = n-1;
			pp1[y * n + x] += 1;	// hist
			qq1[y * n + x] = s[j];	// sig distribution
		// scaled
			x = (p[j] - mn) / mn * n/16 + n/2;
			y = q[j] / mn * n/16 + n/2;
			if (x < 0) x = 0;
			if (x >= n) x = n-1;
			if (y < 0) y = 0;
			if (y >= n) y = n-1;
			pp2[y * n + x] += 1;	// hist
			qq2[y * n + x] = s[j];	// sig distribution
		}
	}
	[sct1 saveAsKOImage:@"IMG_sct1"];	// raw
	[sct2 saveAsKOImage:@"IMG_sct2"];	// scaled
}

// kFiltW
// single 2D image is assumed
RecImage *
kFiltW(RecImage *f)
{
	RecImage		*f1;
    RecLoopControl *lc;
	int				i, j, ii, jj, k, n;
    int             xDim, yDim;
	int				pixSize, dataLength;
	float			*data, *p, *q, x, y;
    float           tmp, w;

    f1 = [f copy];
    lc = [f1 control];
    xDim = [f1 xDim];
    yDim = [f1 yDim];
	pixSize = [f1 pixSize];
	dataLength = [f1 dataLength];

    [lc rewind];
    [lc deactivateXY];
    n = [lc loopLength];
    for (k = 0; k < n; k++) {
        data = [f1 currentDataWithControl:lc];
        for (i = 0; i < yDim; i++) {
            p = data + i*xDim;
            q = p + dataLength;
            // swap hi-lo
            for (j = 0; j < xDim/2; j++) {
                  tmp = p[j]; p[j] = p[j + xDim/2]; p[j + xDim/2] = -tmp;
                  tmp = q[j]; q[j] = q[j + xDim/2]; q[j + xDim/2] = -tmp;
            }
            // filt
            for (j = 0; j < xDim; j++) {
                x = j - xDim/2;
            //    w = cos(x * M_PI / xDim);
                w = exp(-x*x/2000.0);
                p[j] *= w;
                q[j] *= w;
            }
        }
        [lc increment];
    }
	[f1 fft2d:REC_FORWARD];

	return f1;
}


