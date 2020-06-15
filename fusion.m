//
//
//

#import "RecKit.h"
#import "NumKit.h"
#import "RecImageNCI.h"

void		fusion_test();
void		statistic_test();
void		statistic_test_2d();
RecImage	*fusion(RecImage *mg, RecImage *sg, float gain);

int
main()
{
    @autoreleasepool {
		int			mode = 0;	// 0: fusion, 1: statistics 1D, 2: statistics 2D

		switch (mode) {
		case 0 :
			fusion_test();
			break;
	// statistic
		case 1 :
			statistic_test();
			break;
		case 2 :
			statistic_test_2d();
			break;
		}
        return 0;
    }   // autoreleasepool
}

void
fusion_test()
{
	RecImage    *mg, *sig, *col;

	mg = [RecImage imageOfType:RECIMAGE_REAL xDim:256 yDim:256];
	col = [mg colMapImage:1];
	[col saveAsKOImage:@"IMG_col_dipole"];

	mg  = [RecImage imageWithKOImage:@"IMG_img"];
	sig = [RecImage imageWithKOImage:@"IMG_avg"];
	col = fusion(mg, sig, 0.02);
	[col saveAsKOImage:@"IMG_fusion"];
}

float
count_pixel(RecImage *img)
{
	int		i, n, len = [img dataLength];
	float	*p, *q;

	p = [img data];
	q = p + len;

	n = 0;
	for (i = 0; i < len; i++) {
		if (p[i] != 0 || q[i] != 0) n++;
	}
	return (float)n / len;
}

void
statistic_test()
{
	float	*p1, *p2;
	int		i, k;
	int		cnt;
	float	thres = 0.05;	// 5%
	int		navg = 4;
	int		ntry = 128 * 128;
	float	*tval, pval;
	int		*hst;
	int		nbin = 100;

	p1 = (float *)malloc(sizeof(float) * navg);
	p2 = (float *)malloc(sizeof(float) * navg);

	tval = (float *)malloc(sizeof(float) * ntry);
	hst  = (int *)malloc(sizeof(int) * nbin);

// t-test
	cnt = 0;
	for (k = 0; k < ntry; k++) {
		for (i = 0; i < navg; i++) {
			p1[i] = Num_nrml(0.0, 1.0);
			p2[i] = Num_nrml(0.0, 1.0);
		}
		Num_ttest(p1, navg, p2, navg, tval + k, &pval);
		if (pval < thres) cnt++;
	}
	printf("t-test: p < %3.2f : %f (navg/ntry = %d/%d)\n", thres, (float)cnt/ntry, navg, ntry);
// F-test
	cnt = 0;
	for (k = 0; k < ntry; k++) {
		for (i = 0; i < navg; i++) {
			p1[i] = Num_nrml(0.0, 1.0);
			p2[i] = Num_nrml(0.0, 1.0);
		}
		Num_ftest(p1, navg, p2, navg, tval + k, &pval);
		if (pval < thres) cnt++;
	}
	printf("F-test: p < %3.2f : %f (navg/ntry = %d/%d)\n", thres, (float)cnt/ntry, navg, ntry);

//	histogram(hst, nbin, tval, ntry, -3.0, 3.0);
//	for (i = 0; i < nbin; i++) {
//		printf("%d %d\n", i - nbin/2, hst[i]);
//	}

	free(p1);
	free(p2);
	free(hst);
	free(tval);
}

void
statistic_test_2d()
{
	RecImage	*ref, *stm, *pimg;
	RecLoop		*xLp, *yLp, *avgLp, *phsLp;
	float		*p;
	int			nAvg = 4;
	int			i, n;
	int			dim = 128;
	float		pthres = 0.05;
//	float		av, vr;

	system("rm IMG_*");
	xLp = [RecLoop loopWithDataLength:dim];
	yLp = [RecLoop loopWithDataLength:dim];		// ntry = 16000
	avgLp = [RecLoop loopWithDataLength:nAvg];
	phsLp = [RecLoop loopWithDataLength:1];

	ref = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:phsLp, avgLp, yLp, xLp, nil];
	stm = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:phsLp, avgLp, yLp, xLp, nil];
	p = [ref data];
	n = [ref dataLength] * [ref pixSize];

	printf("Nrml (ref)\n");
	for (i = 0; i < n; i++) {
		p[i] = Num_nrml(0.0, 1.0);
	}

	p = [stm data];
	n = [stm dataLength] * [ref pixSize];
	printf("Nrml (stm)\n");
	for (i = 0; i < n; i++) {
		p[i] = Num_nrml(0.0, 1.0);
	}

	[stm saveAsKOImage:@"IMG_stm"];
	[ref saveAsKOImage:@"IMG_ref"];


	printf("t-test\n");
	pimg = [stm pImageWithRef:ref forLoop:avgLp thres:0 ttest:YES];	// p < 0.05
	printf("%3.2f (%3.2f)\n", count_pixel(pimg), pthres);
	[pimg saveAsKOImage:@"IMG_ttest"];

	printf("F-test\n");
	pimg = [stm pImageWithRef:ref forLoop:avgLp thres:0 ttest:NO];	// p < 0.05
	printf("%3.2f (%3.2f)\n", count_pixel(pimg), pthres);
	[pimg saveAsKOImage:@"IMG_ftest"];
/*
4 	0.056946 	0.001221
8	0.057800	0.002075
10	0.059387	0.002808
20	0.065796	0.004822
30	0.066284	0.006104
40	0.064880	0.007935
50	0.067383	0.006592
60	0.065002	0.006531
80	0.066956	0.008972
160	0.065186	0.008179
*/
}

RecImage *
fusion(RecImage *mg, RecImage *sg, float gain)
{
    RecImage    *img;
    int         i, n;
    float       *p, *q, *pp;
    float       *r, *g, *b;
    float       rho, th;

    pp = [mg data];
    n = [mg dataLength];

    p = [sg data];
    q = p + n;

    img = [RecImage imageOfType:RECIMAGE_COLOR withImage:mg];
    // implement saveAsKOImage:, too
    r = [img data];
    g = r + n;
    b = g + n;

    for (i = 0; i < n; i++) {
        rho = p[i] * p[i] + q[i] * q[i];
        rho = sqrt(rho);
        th = atan2(q[i], p[i]) + M_PI * 0.0;

        r[i] = pp[i] + gain * rho * cos(th);
        g[i] = pp[i] + gain * rho * cos(th + M_PI * 2 / 3);
        b[i] = pp[i] + gain * rho * cos(th + M_PI * 4 / 3);
    }

    return img;
}


