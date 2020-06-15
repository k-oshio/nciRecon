//
//	CSF project
//

#import <RecKit/RecKit.h>
#import <NumKit/NumKit.h>
#import "RecImageNCI.h"

NSString	*base = @"/Users/oshio/images/CSF/DWI-rinkan-1";
//NSString	*base = @"/Users/oshio/images/CSF/DWI-rinkan-2";
//NSString	*ser = @"MPG-PE";
NSString	*ser = @"MPG-RO-1";
//NSString	*ser = @"MPG-RO-2";
//NSString	*ser = @"MPG-RO-10AV";
NSString	*path;

void	test1();	// PCA
void	test2();	// phase
void	test3();	// make multi-loop block from 10av set (8-21)
void	test4();	// make multi-loop block from SAG set (11-27)

int
main(int ac, char *av[])
{
//	test1();
//	test2();
//	test3();
	test4();
}

// make PCA lib method
void
test1()
{
	RecImage		*img;
	RecImage		*aimg, *eimg, *uimg, *vimg, *simg;
	Num_mat			*A, *U, *Vt;
	Num_vec			*s;
	int				nr, nc, n;
	RecLoop			*timeLp;
	int				i, tDim, iDim;
	float			*p, *q;
	int				b;

	b = 1000;	// 200, 500, 1000
	path = [NSString stringWithFormat:@"%@/%@/b%d/b%d-%d-mg", base, ser, b, b, b];
	img = [RecImage imageWithKOImage:path];

	timeLp = [img zLoop];
	tDim = [timeLp dataLength];
	iDim = [img xDim] * [img yDim];

// === SVD
	//make A matrix
	aimg = [RecImage imageOfType:RECIMAGE_REAL xDim:iDim yDim:tDim];
	q = [aimg data];
	p = [img data];
	for (i = 0; i < iDim * tDim; i++) {
		q[i] = p[i];
	}
//	[aimg saveAsKOImage:@"IMG_A"];
	nr = [aimg yDim];
	nc = [aimg xDim];
	n = MIN(nr, nc);
	A = Num_im_to_m(aimg);
	// make U, s, Vt
	U = Num_new_mat(nr, n);
	Vt = Num_new_mat(n, nc);
	s = Num_new_vec(n);
	// call LAPACK
//	Num_svd(A, U, s, Vt);	// 5.5 sec #### update arguments

	simg = Num_v_to_im(s);
//	[simg saveAsKOImage:@"IMG_s"];

	uimg = Num_m_to_im(U);
	[uimg trans];
//	[uimg saveAsKOImage:@"IMG_U"];
	path = [NSString stringWithFormat:@"%@/%@/b%d-U", base, ser, b];
	[uimg saveAsKOImage:path];
	
	vimg = Num_m_to_im(Vt);
//	[vimg saveAsKOImage:@"IMG_Vt"];

	// eigen image
	eimg = [RecImage imageOfType:RECIMAGE_REAL xDim:[img xDim] yDim:[img yDim] zDim:tDim];	// -> preserve img dim
	// size should be the same as self
	q = [eimg data];
	p = [vimg data];		// Vt
	for (i = 0; i < iDim * tDim; i++) {
		q[i] = p[i];
	}
//	[eimg saveAsKOImage:@"IMG_E"];
	path = [NSString stringWithFormat:@"%@/%@/b%d-E", base, ser, b];
	[eimg saveAsKOImage:path];
}

void
test2()
{
	int			i, j, b;
	int			btab[3] = {200, 500, 1000};
	RecImage	*img, *phs;
	NSString	*sertab[3] = {@"MPG-RO-1", @"MPG-RO-2", @"MPG-PE"};
	float		scl = M_PI / 10000;

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			b = btab[j];	// 200, 500, 1000
			path = [NSString stringWithFormat:@"%@/%@/b%d/b%d-%d-mg", base, sertab[i], b, b, b];
			img = [RecImage imageWithKOImage:path];
			path = [NSString stringWithFormat:@"%@/%@/b%d/b%d-%d-phs", base, sertab[i], b, b, b];
			phs = [RecImage imageWithKOImage:path];
			[phs multByConst:scl];

			[img makeComplexWithPhs:phs];
			[img pcorr];
			[img checkNeg0];
			path = [NSString stringWithFormat:@"%@/%@/b%d-test-cpx", base, sertab[i], b];
			[img saveAsKOImage:path];
			phs = [img unwrap2d];
			path = [NSString stringWithFormat:@"%@/%@/b%d-test-phs", base, sertab[i], b];
			[phs saveAsKOImage:path];
		}
	}
}

// 10 av axial set
void
test3()
{
	int			b, d;
	RecImage	*img, *phs, *avg;
	float		scl = M_PI / 1800;	//10000;
	NSString	*inSer, *outSer;

// blk : [del, avg, y, x]

	b = 200;	// 200, 500, 1000
//	ser =  @"MPG-RO-10AV";
	inSer = @"MPG-RO-10AV/b200-img";
	outSer = @"MPG-RO-10AV/b200-proc";
	base = @"/Users/oshio/images/CSF/DWI-rinkan-2";

	for (d = 0; d < 7; d++) {
		path = [NSString stringWithFormat:@"%@/%@/b%d-d%d.mg", base, inSer, b, d];
		img = [RecImage imageWithKOImage:path];
		path = [NSString stringWithFormat:@"%@/%@/b%d-d%d.phs", base, inSer, b, d];
		phs = [RecImage imageWithKOImage:path];
		[phs multByConst:scl];
		[img makeComplexWithPhs:phs];
		[img pcorr];
		[img checkNeg0];

	// complex
		path = [NSString stringWithFormat:@"%@/%@/b%d-%d-cpx", base, outSer, b, d];
		[img saveAsKOImage:path];

	// phase
		phs = [img unwrap2d];
		path = [NSString stringWithFormat:@"%@/%@/b%d-%d-phs", base, outSer, b, d];
		[phs saveAsKOImage:path];

	// difference of mag
		[img magnitude];
		avg = [img avgForLoop:[img zLoop]];
		[img subImage:avg];
		path = [NSString stringWithFormat:@"%@/%@/b%d-%d-dif", base, outSer, b, d];
		[img saveAsKOImage:path];
	}
}



