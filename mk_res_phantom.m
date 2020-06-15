//
//  resolution phantom -> move to testProc
//

#import <RecKit/RecKit.h>

RecImage    *ph;

RecImage    *mk_phantom(void);  // stripe
RecImage    *mk_phantom2(void); // stochastic
void		mk_phantom3(void); // complex f-test

int
main()
{
//    ph = mk_phantom();
//    ph = mk_phantom2();
//    [ph saveAsKOImage:@"phantom/PH_res.img"];

	mk_phantom3();

    return 0;
}

RecImage *
mk_phantom(void)
{
    RecImage    *img;
    int         i, j, k, l, n, m;
    int         x0, x, y;
    float       val;
    int         res = 1024;
    int         wd = 800;   // width of non-zero part
    int         hi = 600;   // height of non-zero part
    int         nblk = 10;
    int         blksz = wd / nblk;
    float       *p;

    img = [RecImage imageOfType:RECIMAGE_REAL xDim:res yDim:res];
    p = [img data];

    val = -1;   // background
    for (i = 0; i < hi; i++) {
        y = (res - hi) / 2 + i;
        for (j = 0; j < wd; j++) {
            x = (res - wd) / 2 + j;
            p[y * res + x] = val;
        }
    }
            
    for (k = 0; k < nblk; k++) {   // res block
        x0 = res/2 - wd/2 + k * blksz;   // left edge
        m = k + 1;  // width of single stripe
        for (i = 0; i < hi; i++) {
            n = blksz / m / 2;  // # of stripes
            for (j = 0; j < n; j++) {
                val = 1.0;
                for (l = 0; l < m; l++) {
                    x = x0 + j * m * 2 + l;
                    y = res/2 - hi/2 + i;
                    p[y * res + x] = val;
                }
            }
        }
    }


    return img;
}

RecImage *
mk_phantom2(void)
{
    RecImage    *img;
    int         i, j;
    float       val;
    int         res = 512;
    float       *p;
    RecNoise    *nz = [RecNoise noise];

    img = [RecImage imageOfType:RECIMAGE_REAL xDim:res yDim:res];
    p = [img data];

    for (i = 0; i < res; i++) {
        for (j = 0; j < res; j++) {
            val = [nz nrml];
            p[i * res + j] = val;
        }
    }

    return img;
}

void
mk_phantom3(void) // complex f-test
{
	RecImage	*img;
	int			i, j, k;
	float		val, r;
	int			res = 128;
	int			nAvg = 100;
	int			n;
	float		*p;
    RecNoise    *nz = [RecNoise noise];

    img = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:res yDim:res zDim:nAvg];
    p = [img data];
	n = [img dataLength] * [img pixSize];
    for (i = 0; i < n; i++) {
		val = [nz nrml];
		p[i] = val;
    }
	[img saveAsKOImage:@"phantom/PH_nz1.img"];
    for (i = 0; i < n; i++) {
		val = [nz nrml];
		p[i] = val;
    }
	[img saveAsKOImage:@"phantom/PH_nz2.img"];
    for (k = 0; k < nAvg; k++) {
    for (i = 0; i < res; i++) {
    for (j = 0; j < res; j++) {
		r = (i - res/2)*(i - res/2) + (j - res/2)*(j - res/2);
		r = sqrt(r);
		val = [nz nrml];
		if (r < 10.0) {
			val += 0.4;
		}
		p[k * res * res + i * res + j] = val;
    }
	}
	}
	[img saveAsKOImage:@"phantom/PH_nz3.img"];
}


