//
//
//

#include <stdio.h>
#include <math.h>

void test0();   // ??
void test1();   // vector filter test
void test2();   // muliti-angle
void test3();   // BOLD

int
main()
{
//    test0();    // ??
//    test1();
//    test2();
    test3();

    return 0;
}

// ???
void
test0(void)
{
	int		i, ix;
	float	m, x, y, th;
	char	buf[256];

	for (i = 0; i < 360; i++) {
		if (gets(buf) == NULL) break;
		sscanf(buf, "%d %f", &ix, &m);
		th = i * M_PI / 180;
	//	m -= 2000.0;
		x = m * cos(th);
		y = m * sin(th);
		printf("%f %f\n", x, y);
	}
}

#define N 90

void
test1()
{
    int     i;
    float   th, c, s;
    float   A, B, C, D;

    for (i = 0; i <= N; i++) {
        th = i * 2.0 * M_PI / N;
        c = cos(th);
        s = sin(th);
        A = c * c;
        B = s * c;
        C = 1/sqrt(2) * (s + c) * c;
        D = 1/sqrt(2) * (s - c) * c;
        printf("%f %f %f %f %f\n", th * 180 / M_PI, A, B, C, C - D);
    }
        
}
void
test2()
{
    int     i, j;
    float   th, th2;
    float   A, B, S, C, D;
    float   sum_c, sum_s;

    for (i = 0; i <= N; i++) {
        th = i * 2.0 * M_PI / N;
        A = cos(th);
        B = sin(th);
        sum_c = sum_s = 0;
        for (j = 0; j < 16; j++) {
            th2 = th * j;
            if (j % 2 == 0) {
                C = 1;
                D = 0;
            } else {
                C = 1;
                D = 0;
            }
            sum_c += A * (C * cos(th2) + D * sin(th2));
            sum_s += B * (C * cos(th2) + D * sin(th2));
        }
//        printf("%f %f %f %f %f\n", th * 180 / M_PI, A, B, sum_c, sum_s);
        printf("%f %f %f\n", th * 180 / M_PI, sum_c, sum_s);
    }
        
}

void
test3()
{
    int     i, j;
    float   f;
    float   th, ph, r;
    float   r2, ph2;
    float   C1, C2, S1, S2;
    float   cp, ct, sp, st;
    float   A = 1.0;

    for (i = 0; i <= N; i++) {
        ph = i * 2.0 * M_PI / N;
        printf("%f ", ph);
        
        cp = cos(ph)*cos(ph);
        sp = sin(ph)*sin(ph);

        for (j = 0; j < 10; j++) {
            th = 0.05 * M_PI * j;

            ct = cos(th)*cos(th);
            st = sin(th)*sin(th);

            S1 = st*st;
            C1 = cp * ct + sp;
            C2 = cp * ct - sp;
            f = A * S1 * C2 / (C1 * C1);
            printf("%f ", f);
        }
        printf("\n");
    }
        
}



