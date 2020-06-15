//
//	dipole imaging
//  with k-space inverse filtering
//

#import <RecKit/RecKit.h>


RecImage		*kFilt(RecImage *f, float dir);	// dipole model
void            toCsv(RecImage *f, RecImage *mask);
RecImage        *toImg(char *fname, RecImage *mask);
void            plotCsv(char *csvPath, char *datPath);

int
main(int ac, char *av[])
{
	RecImage	*f, *fx, *fy, *mag, *mask, *mn;
    RecImage    *ff[256], *tmp;
    int         nang;
    RecLoop     *zLp;
	char		fname[256];
	NSString	*path;
	float		*p, *q, re, im, r2;
	int			i, len;
	int			nphs = 2, sliceNum = 0, slicePos = 0;
	int			nslc = 5;

    if (ac <= 1) {
        printf("kFilt <sliceNum>\n");
    //    exit(0);
    }
	if (ac > 1) {
        sliceNum = atoi(av[1]);
        slicePos = sliceNum * 2;
        if (slicePos >= nslc) {
            slicePos -= nslc;
        }

		path = [NSString stringWithFormat:@"raw%03d.block", sliceNum];
		f = [RecImage imageWithKOImage:path];
	} else {
	//	f = [RecImage imageWithKOImage:@"raw0818.block"];
	//	nphs = 2;
		f = [RecImage imageWithKOImage:@"PH_k.img"];
		nphs = 0;
	}
	if (f == nil) {
		printf("Couldn't open image file [%s]\n", fname);
		exit(0);
	}
	if (ac > 2) {
		nphs = atoi(av[2]);
	}
	mag = [f copy];
	[mag fft2d:REC_FORWARD];
    mask = [mag sumForLoop:[mag topLoop]];
    [mask thresAt:0.05];
	[mask saveAsKOImage:@"IMG_mask"];
    
    mag = [mag sumForLoop:[mag topLoop]];
	mag = [mag oversample];
    path = [NSString stringWithFormat:@"IMG_mag%d", slicePos];
	[mag saveAsKOImage:path];

// calc inverse for multiple directions
    nang = 8;
    for (i = 0; i < nang; i++) {
        ff[i] = kFilt(f, i * M_PI/nang);	// dipole model (th =   0, x)
    }
// combine and make x/y
    // 45
    fx = [RecImage imageWithImage:ff[0]];
    for (i = 0; i < nang; i++) {
        r2 = cos(M_PI/4 + i * M_PI/nang);
        tmp = [ff[i] copy];
        [tmp multByConst:r2];
        [fx addImage:tmp];
    }

    // 135
    fy = [RecImage imageWithImage:ff[0]];
    for (i = 0; i < nang; i++) {
        r2 = cos(3*M_PI/4 + i * M_PI/nang);
        tmp = [ff[i] copy];
        [tmp multByConst:r2];
        [fy addImage:tmp];
    }

// make single image
	q = [fx data];
	p = [fy data];
	len = [fx dataLength];
	q += len;
	for (i = 0; i < len; i++) {
		q[i] = p[i];
	}
    if (nphs == 0) {
        [fx saveAsKOImage:@"IMG_PH_xy"];
        return 0;    // for numerical phantom
    }

    mn = [fx copy];
    zLp = [mn topLoop];
	mn = [mn avgForLoop:zLp];
	[mn addLoop:zLp];
	[fx subImage:mn];

    toCsv(fx, mask);
	fx = [fx oversample];
    path = [NSString stringWithFormat:@"IMG_xy%d", slicePos];
	[fx saveAsKOImage:path];

    printf("kFilt ok\n");

    system("echo calling fastICA\n");
    system("/usr/local/bin/octave /Users/oshio/Projects6/RecKit/nciRecon/test_ica.m");
    // csv to image
    fx = toImg("IMG_s.csv", mask);
	fx = [fx oversample];
    path = [NSString stringWithFormat:@"IMG_ica%d", slicePos];
	[fx saveAsKOImage:path];
    // csv to plot file
    sprintf(fname, "ICA_a%d.dat", slicePos);
    plotCsv("ICA_a.csv", fname);

	return 0;
}

RecImage *
kFilt(RecImage *f, float dir)
{
	RecImage		*f1, *f2;
	RecImage		*mag;
	RecLoopControl	*lc;
	int				xDim = [f xDim];
	int				yDim = [f yDim];
	int				i, j, k, ix;
	int				n, dataLength;
	float			*p, *q, *pp, mg;
	float			x, y, r, cs, tth, kr, w, invr;
	int				mode = 4;       // 0: cos, 1: 1/cos 2: 1/(cos + s) 3: r / cos 4: 3-d
	float			mgs = 0.01;     //                              0.001
	float			sig = 0.30;		// noise in imag part [0..1]    0.20
	float			sd = 30.0;		// high pass freq               30.0
    float           lsd = 80.0;     // low pass freq                80.0

	dataLength = [f dataLength];
	lc = [f control];
	[lc deactivateXY];
	n = [lc loopLength];

// mag
	mag = [f copy];
	[mag fft2d:REC_FORWARD];
	[mag magnitude];
	mgs *= [mag maxVal];

// take imag part
	f1 = [f copy];
	[f1 fft2d:REC_FORWARD];
	p = [f1 data];
	for (i = 0; i < dataLength; i++) {
		p[i] = 0;
	}
	[f1 fft2d:REC_INVERSE];

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
                tth = atan2(y, x);
                cs = cos(tth + dir);
                
                // div type
				switch (mode) {
				case 0:
 					invr = cs;
                    w = 1.0 - exp(-r*r/(sd*sd));	// -> inv gauss hi-pass
					p[ix] *= w * invr;
					q[ix] *= w * invr;
					break;
				case 1:
					if (cs > sig) {
						invr = 1/cs;
					} else {
						invr = 0;
					}
                    w = 1.0 - exp(-r*r/(sd*sd));	// -> inv gauss hi-pass
					p[ix] *= w * invr;
					q[ix] *= w * invr;
					break;
				case 2:
					invr = cs / (cs*cs + sig*sig);
                    w = 1.0 - exp(-r*r/(sd*sd));	// -> inv gauss hi-pass
					p[ix] *= w * invr;
					q[ix] *= w * invr;
					break;
				case 3:
					invr = r * cs / (cs*cs + sig*sig);
                    w = exp(-r*r/(lsd*lsd));	// -> inv gauss hi-pass
					p[ix] *= w * invr;
					q[ix] *= w * invr;
					break;
				case 4: // 3-d psf
                    kr = 1.0 / (0.5 * exp(-r / 10.0) + 0.5 * exp(-r / 50.0));
					invr = kr * cs / (cs*cs + sig*sig);
               //     if (r >= 96) {
               //         w = 0;
               //     } else {
               //         w = cos(r * M_PI / 2 / 90.0);  // half sin
               //     }
                    w = 1.0;
 					p[ix] *= w * invr;
					q[ix] *= w * invr;
					break;
				}
			}
		}
		[lc increment];
	}

// FT
	[f1 fft2d:REC_FORWARD];

/* take phase */
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
	return f1;
}

// convert to csv format for Octave processing (ICA)
void
toCsv(RecImage *f, RecImage *mask)
{
    int     i, j, ix;
    int     imgSize = [mask dataLength];
    int     nValidPix;
    int     nSlice = [f nImages];
    int     dataLength = [f dataLength];
    float   *p, *q, *msk;
    FILE    *fp;

    fp = fopen("IMG_x.csv", "w");

    p = [f data];
    q = p + [f dataLength];
    msk = [mask data];
    nValidPix = 0;
    for (i = 0; i < imgSize; i++) {
        if (msk[i] != 0) nValidPix++;
    }

    for (i = 0; i < nSlice; i++) {
        for (j = 0; j < imgSize; j++) {
            if (msk[j] != 0) {
                fprintf(fp, "%f,", p[j]);
            }
        }
        for (j = ix = 0; j < imgSize ; j++) {
            if (msk[j] != 0) {
                if (ix >= nValidPix - 1) {
                    fprintf(fp, "%f", q[j]);
                    ix++;
                } else {
                    fprintf(fp, "%f,", q[j]);
                    ix++;
                }
            }
        }
        if (i < nSlice - 1) {
            fprintf(fp, "\n");
        }
        p += imgSize;
        q = p + dataLength;
    }
    fclose(fp);
}

RecImage *
toImg(char *fname, RecImage *mask)
{
	RecImage	*f;
	RecLoop		*kx, *ky, *kz;
	float		*p, *q, *msk;
	int			i, j, ix, c;
    int         imgSize, nSlice, nEl, dataLength;
    FILE        *fp;

    fp = fopen(fname, "r");   // output (ICA)
//    fp = fopen("IMG_x.csv", "r"); // input (for testing)
    // count # of lines
    nSlice = 0;
    while ((c = getc(fp)) != EOF) {
        if (c == '\n') nSlice++;
    }

    msk = [mask data];
    dataLength = [mask dataLength];

    nEl = 0;
    for (i = 0; i < dataLength; i++) {
        if (msk[i] != 0) nEl++;
    }
    printf("nSlice = %d, nEl = %d\n", nSlice, nEl);

	kx = [RecLoop loopWithName:@"kx" dataLength:128];
	ky = [RecLoop loopWithName:@"ky" dataLength:128];
	kz = [RecLoop loopWithName:@"kz" dataLength:nSlice];
	f = [RecImage imageOfType:RECIMAGE_COMPLEX withLoops:kz, ky, kx, nil];
    dataLength = [f dataLength];
    imgSize = [mask dataLength];


    p = [f data];
    q = p + [f dataLength];
    msk = [mask data];

    rewind(fp);
    for (i = 0; i < nSlice; i++) {
        for (j = ix = 0; j < imgSize; j++) {
            if (msk[j] != 0) {
                fscanf(fp, "%f,", &p[j]);
                ix++;
            }
            if (ix >= nEl) break;
        }
        for (j = ix = 0; j < imgSize ; j++) {
            if (msk[j] != 0) {
                fscanf(fp, "%f,", &q[j]);
                ix++;
            }
           if (ix >= nEl) break;
        }
        fscanf(fp, "\n");
        p += imgSize;
        q = p + dataLength;
    }
    fclose(fp);

    return f;
}

void
plotCsv(char *csvPath, char *datPath)
{
    FILE    *fp;
	int		i, j, n, ix;
    int     xdim, ydim;
    int     c;
    float   *buf, *p, val;
    float   len;

    fp = fopen(csvPath, "r");
    if (fp == NULL) return;

    xdim = 0;
    while (((c = getc(fp)) != EOF) && (c != '\n')) {
        if (c == ',') xdim++;
    }
    xdim += 1;

    rewind(fp);
    ydim = 0;
    while ((c = getc(fp)) != EOF) {
        if (c == '\n') ydim++;
    }
    
    printf("xdim = %d, ydim = %d\n", xdim, ydim);

// read csv
    buf = (float *)malloc(xdim * ydim * sizeof(float));
    rewind(fp);
    p = buf;
    for (i = 0; i < ydim; i++) {
        for (j = 0; j < xdim; j++) {
            fscanf(fp, "%f,", p + j);
         }
         fscanf(fp, "\n");
         p += xdim;
    }
    fclose(fp);

// write .dat
    fp = fopen(datPath, "w");
    for (i = 0; i < ydim; i++) {
        fprintf(fp, "%d ", i);
        for (j = 0; j < xdim; j++) {
            val = *(buf + i * xdim + j);
            fprintf(fp, "%7.4f ", val);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    free(buf);
}

