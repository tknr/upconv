
#ifndef FFT_FILTER_H
#define FFT_FILTER_H

#ifndef M_PI
#define M_PI		(3.14159265358979323846)
#endif

#ifndef DWORD
#define DWORD		unsigned long
#endif

#ifndef TRUE
#define TRUE		(1)
#endif

#ifndef FALSE
#define FALSE		(0)
#endif

// サンプルを処理するデータ型
#define SSIZE	signed long long int

typedef struct {
	int err;
	long inSampleR;
	long outSampleR;
	long hfc;
	long lfc;
	int  analyze_mode;
	int  disable_eq;
	int  hi_sample_mode;
	int  eq_flag;
	int  lpf_flag;
	int  lfa_flag;
	int  cut_high_dither;
	int  lvadj_flag;
	double *eq_ref_max;
	double *eq_ref_avg;
	double *eq;
	double *lfa_eq;
	long   eq_ref_count;
	int deEmphasis;
	int dsd_fmt;
	char *abort_filename;
	int abort_percent;
} FFT_PARAM;

void fftFilter(DWORD inSample,DWORD outSample,FIO *fp,FIO *tmp,FFT_PARAM *param);
void fftFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,FFT_PARAM *param,int id);
void copyToFFTW(fftw_complex *fftw,SSIZE *buf,long size);
void windowFFTW(fftw_complex *fftw,long size);
void cutFFTW(fftw_complex *fftw,long index,long size);
void *al_malloc(SSIZE size);
void *al_free(void *ptr);

#endif
