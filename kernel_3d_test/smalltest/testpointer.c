#include <stdio.h>

int main()
{


	float * a = (float*)malloc(sizeof(float)*10);


	float * b = (float*)malloc(sizeof(float)*10);

	float *tmp = NULL;
	int i;
	for(i=0;i<10;i++){
		a[i] = i;
		b[i] = i+100;
	}

	
	for(i=0;i<10;i++){
		fprintf(stderr, "%f\t", a[i]);
	}
	fprintf(stderr, "\n");
	for(i=0;i<10;i++){
		fprintf(stderr, "%f\t", b[i]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");
	
	tmp = a;
	a = b;
	b = tmp;

	for(i=0;i<10;i++){
		fprintf(stderr, "%f\t", a[i]);
	}
	fprintf(stderr, "\n");
	for(i=0;i<10;i++){
		fprintf(stderr, "%f\t", b[i]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");
	
	return 1;
}
