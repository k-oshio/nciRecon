//
//

#import <RecKit/RecKit.h>

int
main(int ac, char *av[])
{
	NSArray			*lines;
    NSString		*content, *path;
	char			cwd[256];
    NSError			*err;
	int				i, j, n, siz;
	float			val, *delay;
	int				d1, d2;
	RecImage		*img;
	RecLoopControl	*lc;
	float			*p;
	FILE			*fp;

//    content = [NSString stringWithContentsOfFile:@"meas_files/Chckerboard_event_25_75ms_Cum-1019-1.xlsx.csv"
//        encoding:NSASCIIStringEncoding error:&err];
	if (ac < 2) {
		printf("readCsv <csvfile>\n");
		exit(0);
	}
	path = [NSString stringWithFormat:@"%s", av[1]];
	if (getcwd(cwd, 256) != NULL) {
		path = [NSString stringWithFormat:@"%s/%@", cwd, path];
	}

    content = [NSString stringWithContentsOfFile:path encoding:NSASCIIStringEncoding error:&err];

    lines = [content componentsSeparatedByString:@"\r\n\r\n"];	// windows (CRLF CRLF)
	n = (int)[lines count];
	printf("%d lines found\n", n);

	delay = (float *)malloc(sizeof(float) * n);

	d1 = d2 = 0;
	for (i = 0; i < n; i++) {
		val = delay[i] = [[lines objectAtIndex:i] floatValue];
		if (val < 50) {
			d1++;
		} else {
			d2++;
		}
	}
	printf("===========\n");
	printf("delay1 : %d, delay2 : %d\n", d1, d2);

//	make test image
	n = 10 + n * 2;	// dda + off/delay
	img = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:128 yDim:128 zDim:n];
	lc = [img control];
	siz = [img xDim] * [img yDim];

	[lc rewind];
	[lc deactivateXY];

	fp = fopen("protocol.dat", "w");
	for (i = 0; i < n; i++) {
	//
//printf("zpos = %d\n", [lc zPosition]);
		if (i < 10) {
			fprintf(fp, "-50\n");
		} else
		if (i % 2 == 0) {
			// off
			p = [img currentDataWithControl:lc];
			for (j = 0; j < siz; j++) {
				p[j] = 50;
			}
			fprintf(fp, "0\n");
		} else {
			// delay (random)
			p = [img currentDataWithControl:lc];
			val = delay[(i - 10) / 2];
			for (j = 0; j < siz; j++) {
				p[j] = val;
			}
			fprintf(fp, "%5.2f\n", val);
		}
		[lc increment];
	}
	[img saveAsKOImage:@"test_rand.img"];	// IMGxxx is removed at each iteration

	free(delay);
}