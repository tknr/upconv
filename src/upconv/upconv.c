//---------------------------------------------------------------------------
/****************************************************************************/
/* upconv(C) 2007-2019 By 59414d41											*/
/*																			*/
/*																			*/
/****************************************************************************/

/*--- Log ------------------------------------------------------------------
 * Ver 0.10 <07/03/15> - 作成
 * Ver 0.20 <07/04/19> - DFT バージョン追加
 * Ver 0.30 <08/08/17> - 新方式に変更
 * Ver 0.40 <09/02/20> - fftw バージョンに変更しいったんFix
 * Ver 0.50 <09/05/18> - 高域の補間方法をいったんFix
 * Ver 0.60 <09/06/02> - BCB6とVC6に対応
 * Ver 0.61 <09/06/06> - メモリマップドファイルを使用しないようにした
 * Ver 0.70 <09/06/28> - 処理方法を変更
 * Ver 0.80 <09/07/10> - hfa2 の処理を改良しマルチコアに対応するために処理を分離
 * Ver 0.90 <09/09/24> - メモリ使用量削減
 * Ver 0.91 <09/09/27> - ダウンサンプリング時のバグ修正、ノーマライズ時のバグ修正
 * Ver 0.92 <09/09/29> - hfa2の補間方法を変更
 * Ver 0.93 <09/10/06> - hfa2の補間処理を変更
 * Ver 0.94 <09/11/01> - hfa3追加、バグ修正、パラメータファイルの採用
 * Ver 0.95 <09/11/11> - hfa3にノイズとのブレンド具合を指定可能にした
 * Ver 0.96 <09/11/15> - hfa補間時の周波数特性を修正
 * Ver 0.97 <09/11/22> - ビット拡張処理追加、hfa補間時の周波数特性を指定可能にした
 * Ver 0.98 <10/01/11> - デエンファシスに対応
 * Ver 0.99 <10/03/14> - hfa3の補間方法を変更
 * Ver 1.00 <10/03/14> - GCC(MinGW)に対応
 * Ver 1.01 <10/04/11> - OpenMPに対応
 * Ver 1.02 <10/07/26> - スピーカーごとの調整機能追加
 * Ver 1.03 <10/09/14> - hfc autoに対応
 * Ver 1.04 <10/11/02> - スピーカーごとの調整機能バグ修正
 * Ver 1.05 <10/12/27> - イコライザ機能修正
 * Ver 1.06 <11/01/07> - lfa 対応、ソースコードの整理
 * Ver 1.07 <11/10/01> - テンポラリファイル対応
 * Ver 1.08 <12/02/28> - fio 対応
 * Ver 1.09 <13/04/27> - hfa3の補間方法変更
 * Ver 1.10 <18/10/06> - サンプリングレート追加
 * Ver 1.20 <19/10/12> - いろいろな改良
 */

// [仕様]
// upconv.exe 入力ファイル名 出力ファイル名 デフォルトパラメータファイル名 パラメーター
//
// [作業用ファイル]
// ◇ 出力ファイル名.param
// 1行目は元のwavファイル or 元のmp3ファイル or flacファイル or wavpackファイルの名前
// ノーマライズ用音量データ(Ch1)
// ノーマライズ用音量データ(Ch2)
// ...

// ◇ 出力ファイル名.files
// 作業用に作成したファイル名
// ...

// ◇ 出力ファイル名.r1
// 変換中の1Ch分のデータ(64bit,有効60bitのrawデータ)
// ...

// ◇出力ファイル名.r1.param
// 1Ch分の変換後データのノーマライズ用音量データ

#define STR_COPYRIGHT	"upconv.exe (c) 2019 Ver 1.20 By 59414d41\n\n"
#define STR_USAGE		"upconv.exe in-file out-file def_paramfile parameter\n"

#if 0
#define	PRINT_LOG(s)	do {																	\
							FILE *log;															\
							log = fopen("d:\\upconv.log","a");									\
							if (log) {															\
								fprintf(log,"%s [%d] %s\n",__FUNCTION__,__LINE__,s);			\
								fclose(log);													\
							}																	\
						} while (0)
#else
#define	PRINT_LOG(s)	//
#endif


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <unistd.h>
#include "upconv.h"
#include "fileio.h"
#include "fftw3.h"
#include "fft_filter.h"
#include "./../PLG_AUDIO_IO/PLG_AudioIO.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// 最大サンプリング数
#define MAX_SAMPLE_N	(1536000*2)

// ノーマライズ情報
typedef struct {
	SSIZE	min;
	SSIZE	max;
	SSIZE	avg;
} NORM_INFO;

typedef struct {
	long cutOff[257];
} BITEXT_INFO;

typedef struct {
	int upconv;
	int err;
	int errLine;
	int disRc;
	int disEq;
	long inSampleR;
	long outSampleR;
	long validOutSampleR;
	int hfa;
	int low_adjust;
	int overSamp;
	int abe;
	int abeNX;
	int abeFnLevel;
	int cutLowData;
	int smLowData;
	int post_abe;
	int adaptiveFilder;
	int w;
	int iw;
	int channel_count;
	int hiDither;
	int hfc_auto;
	int cut_high_dither;
	int enable_hfc;
	long hfc;
	int enable_lfc;
	long lfc;
	long lpf;
	int enable_nr;
	long nr;
	long fio;
	int nrLV;
	int Adjust_enable;
	// hfa2,3 option
	int hfa_preset;
	int sig1Enb,sig2Enb,sig3Enb;
	int sig1AvgLineNx;
	int sig1Phase,sig2Phase,sig3Phase;
	int hfaNB;
	int hfaFast;
	int hfaWide;
	int hfaDiffMin;
	// 
	int deEmphasis;
	BITEXT_INFO beInfo;
	char tempPath[4];
	int thread;
	int downSample;
	int ana;
	int sp_ana;
//	double *ana_pw;
//	double *ana_avgpw;
//	double *eq_pw;
	int    eq_pwcount;
	long   ana_avgcount;
//	double *sp_eq;
//	double *eq;
	char *sp_path;
	char *nd_path;
	int r1_flag;
	int eq_flag;
	int enable_filter;
	double hfa3_max;
	int pwAdj;
	int dsd_fmt;
	int mc_flag;
	char abort_filename[_MAX_PATH];
	int abort_percent;
	int cpu_pri;
} PARAM_INFO;

typedef struct {
	double	power[65536*2];
	double	phase[65536*2];
	double	pw[65536*2];
	double	diff[65536*2];
	double	avg[65536*2];
	double	base[65536*2];
	double	baseToAdj[65536*2];
	int		sign[65536*2];
	int		pw_cnt[65536*2];
	int		nSample;
	int		validSamplingRate;
	long	samplingRate;
	int		do2Idx[360];
	double	*phaseX;
	double	*phaseY;
	int log;
} OVERTONE_INFO;

SSIZE *diskBuffer;
SSIZE *diskBuffer2;
NORM_INFO NormInfo;

/*--- wav2raw Function Prototype ---------------------------------------------------*/
extern int to_raw_main(int argc, char *argv[]);
extern int dsf_main(int argc, char *argv[]);
extern int mc_main(int argc,char *argv[]);
extern int to_wav_main(int argc, char *argv[]);
extern int start_exec(int argc,char *argv[],int cpu_pri,HANDLE *ret_hStdOutRead,HANDLE *ret_process_id,HANDLE *ret_thread_id);
//extern pid_t fork(void);

/*--- Function Prototype ---------------------------------------------------*/
void SamplingRateConvert(char *rawFile,PARAM_INFO *param);
void anaLFA_Param(FFT_PARAM *param);
void anaHFC_AutoParam(FFT_PARAM *param);
void spAnalyze(SSIZE inSample,FIO *fp,PARAM_INFO *param);
void adjBitExtension(SSIZE inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param);
void genNoise(long hfc,SSIZE inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param);
void genOverTone(long hfc,SSIZE inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param,FFT_PARAM *fft_param);
void genOverToneSub(long hfc,SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,OVERTONE_INFO *ovInfo,PARAM_INFO *param,FFT_PARAM *fft_param);
//int bpFilter(long lfc,long hfc,DWORD inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param);
//void bpFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,long lfc,long hfc,PARAM_INFO *param);

void noiseCut(long nfs,SSIZE inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param);
void anaOverToneHFA2(OVERTONE_INFO *ovInfo,PARAM_INFO *param);
void anaOverToneHFA3(OVERTONE_INFO *ovInfo,PARAM_INFO *param);
void outTempFile(FIO *fio,void *in,SSIZE size,PARAM_INFO *param);
double normalNoise(void);
void adjPinkFilter(int mode,long fftSizeOut,fftw_complex *fftw_out2,PARAM_INFO *param);
void merageTempFile(char type,int normFlag,FIO *in1,FIO *in2,FIO *out,SSIZE inSample,PARAM_INFO *param);
void copyToFFTW(fftw_complex *fftw,SSIZE *buf,long size);
void windowFFTW(fftw_complex *fftw,long size);
void cutFFTW(fftw_complex *fftw,long index,long size);
static int chkAbort(PARAM_INFO *param,int percent,int diff);

//---------------------------------------------------------------------------
// Function   : upconv_main
// Description: upconv メイン関数
//
//
int main(int argc, char *argv[])
{
	FILE *fp;
	FILE *fp_files;
	char workpath[_MAX_PATH];
	char workstr[2048];
	char tmppath[_MAX_PATH];
	char workdrive[_MAX_DRIVE];
	char workdir[_MAX_DIR];
	char workfname[_MAX_FNAME];
	char workext[_MAX_EXT];
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];
	char work[2048];
	char pparam[2048];
	char ddBuf[100];
	long paramFlag;
	char *rn_fname;
	char *p1,*p2;
	double dd;
	long i;
	PARAM_INFO param;
	STARTEXEC_INFO startexec_info[6];
	long temp,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9;
	int retCode = 0;
	SSIZE max,avg;
	double percent;
	int c;
	int retry;

	do {
		memset(&param,0,sizeof (PARAM_INFO));
		memset(startexec_info,0,sizeof (STARTEXEC_INFO) * 6);

		diskBuffer	= (SSIZE *)calloc(4 * 1024 * 1024L,sizeof (SSIZE));
		diskBuffer2 = (SSIZE *)calloc(4 * 1024 * 1024L,sizeof (SSIZE));
		if (diskBuffer == NULL || diskBuffer2 == NULL) {
			param.err = STATUS_MEM_ALLOC_ERR;param.errLine = __LINE__;
			break;
		}
		memset(&NormInfo,0,sizeof (NORM_INFO));

		param.ana_avgcount = 0;
		param.err = STATUS_SUCCESS;
		param.disRc = 0;
		param.disEq = 0;
		param.enable_hfc = 0;
		param.hfc = -1;
		param.enable_lfc = 0;
		param.lfc = -1;
		param.lpf = -1;
		param.enable_nr = 0;
		param.nr = -1;
		param.nrLV = 0;
		param.hfaNB = 0;
		param.thread = 1;
		param.sp_ana = 0;
		param.hiDither = 0;
		param.hfc_auto = 0;
		param.r1_flag = 0;
		param.eq_flag = 0;
		param.post_abe = 0;
		param.fio = -1;
		param.abeNX = 1;
		param.adaptiveFilder = 1;
		param.pwAdj = 0;
		param.cut_high_dither = 0;
		param.enable_filter = 0;
		param.Adjust_enable = 0;
		param.cutLowData = 0;
		param.overSamp = 0;
		param.dsd_fmt  = -1;
		param.hfa_preset = 1;
		param.sig1Enb = 1;
		param.sig2Enb = 1;
		param.sig3Enb = 1;
		param.sig1AvgLineNx = 3;
		param.sig1Phase = -4;
		param.sig2Phase = -1;
		param.sig3Phase = -3;
		param.hfaNB = 0;
		param.upconv   = 0;	// upconv 使用フラグ
		param.channel_count = 0;
		param.mc_flag = 0;
		param.abort_percent = 0;
		param.cpu_pri = 0;
		paramFlag = 0;
		pparam[0] = '\0';

		sprintf(work,"argc:%d,argv:[%s,%s,%s,%s,%s]",argc,argv[0],argv[1],argv[2],argv[3],argv[4]);
		PRINT_LOG(work);
		if (argc == 5) {
			// default parameter ファイル
			fp = fopen(argv[3],"r");
			if (fp == NULL) {
				retCode = STATUS_PARAMETER_ERR;param.errLine = __LINE__;
				break;
			}
			
			// パラメータの読みこみ
			if (fgets(work,2047,fp) == NULL) {
				retCode = STATUS_PARAMETER_ERR;param.errLine = __LINE__;
				break;
			}
			p1 = strrchr(work,'\n');if (p1 != NULL) *p1 = '\0';
			strcat(pparam,work);
			strcat(pparam," ");
			if (strlen(argv[4]) >= 1) strcat(pparam,argv[4]);

			// tmpファイル用の作業ディレクトリ
			if (fgets(workpath,_MAX_PATH,fp) == NULL) {
				retCode = STATUS_PARAMETER_ERR;param.errLine = __LINE__;
				break;
			}
			p1 = strrchr(workpath,'\n');if (p1 != NULL) *p1 = '\0';
			if (strlen(workpath) >= 2 && workpath[strlen(workpath) - 1] != '\\') strcat(workpath,"\\");

			fclose(fp);
			fp = NULL;

			p1 = pparam;
			p2 = strchr(p1,(int)' ');
			for (;p1 != NULL;) {
				if (p2 != NULL) {
					*p2 = '\0';
				}

				if (sscanf(p1,"-s:%ld",&temp) == 1) {
					switch (temp) {
						case 32000:
						case 44100:
						case 48000:
						case 88200:
						case 96000:
						case 176400:
						case 192000:
						case 352800:
						case 384000:
						case 705600:
						case 768000:
							paramFlag |= 0x0002;
							param.outSampleR = temp;
							param.validOutSampleR = temp;
							break;
					}
				}
				if (sscanf(p1,"-is:%ld",&temp) == 1) {
					switch (temp) {
						case 22050:
						case 24000:
						case 32000:
						case 44100:
						case 48000:
						case 88200:
						case 96000:
						case 176400:
						case 192000:
						case 352800:
						case 384000:
						case 705600:
						case 768000:
							param.inSampleR = temp;
							paramFlag |= 0x0001;
							break;
					}
				}
				if (sscanf(p1,"-ms:%ld",&temp) == 1) {
					switch (temp) {
						case 32000:
						case 44100:
						case 48000:
						case 88200:
						case 96000:
						case 176400:
						case 192000:
						case 352800:
						case 384000:
						case 705600:
						case 768000:
							paramFlag |= 0x0002;
							param.outSampleR = temp;
							param.validOutSampleR = temp;
							break;
					}
				}
				if (sscanf(p1,"-w:%ld",&temp) == 1) {
					if (temp == 16 || temp == 24 || temp == 32) {
						param.w = (int)temp;
					}
				}
				if (sscanf(p1,"-iw:%ld",&temp) == 1) {
					if (temp == 16 || temp == 24 || temp == 32) {
						param.iw = (int)temp;
					}
				}
				if (sscanf(p1,"-hfa:%ld",&temp) == 1) {
					switch (temp) {
						case 0:
						case 1:
						case 2:
						case 3:
							paramFlag |= 0x0004;
							param.hfa = (int)temp;
							break;
					}
				}
				if (sscanf(p1,"-enable_hfc:%ld",&temp) == 1) {
					if (temp == 1) param.enable_hfc = 1;
				}
				if (sscanf(p1,"-hfc:%ld",&temp) == 1) {
					if (temp >= 1000 && temp <= (384000 / 2)) {
						param.hfc = temp;
					}
				}
				if (sscanf(p1,"-enable_lfc:%ld",&temp) == 1) {
					if (temp == 1) param.enable_lfc = 1;
				}
				if (sscanf(p1,"-lfc:%ld",&temp) == 1) {
					if (temp >= 0 && temp <= (384000 / 2)) {
						param.lfc = temp;
					}
				}
				if (strcmp(p1,"-abe") == 0) {
					param.abe = 1;
				}
				if (strcmp(p1,"-post_abe") == 0) {
					param.post_abe = 1;
				}
				if (sscanf(p1,"-abe_option:%d,%d,%d",&temp,&temp2,&temp3) == 3) {
					if (temp == 1) param.cutLowData = 1;
					if (temp2 == 1) param.smLowData = 1;
					if (temp3 == 1) {
						param.abeFnLevel = 2;
					} else if (temp3 == 2) {
						param.abeFnLevel = 5;
					} else if (temp3 == 3) {
						param.abeFnLevel = 10;
					} else if (temp3 == 4) {
						param.abeFnLevel = 15;
					} else if (temp3 == 5) {
						param.abeFnLevel = 20;
					}
				}

				if (strcmp(p1,"-cut_high_dither") == 0) {
					param.cut_high_dither = 1;
				}
				if (strcmp(p1,"-low_adjust") == 0) {
					param.low_adjust = 1;
				}
				if (strcmp(p1,"-hfc_auto") == 0) {
					param.hfc_auto = 1;
				}
				if (sscanf(p1,"-oversamp:%ld",&temp) == 1) {
					if (temp >= 0 || temp <= 3) {
						param.overSamp = (int)temp;
					}
				}

#if 0
				if (sscanf(p1,"-lpf:%ld",&temp) == 1) {
					if (temp >= 1000 && temp <= (384000 / 2)) {
						param.lpf = temp;
					}
				}
#endif
				if (sscanf(p1,"-enable_nr:%d",&temp) == 1) {
					if (temp == 1) param.enable_nr = 1;
				}
				if (sscanf(p1,"-nr:%ld",&temp) == 1) {
					if (temp >= 100 && temp <= (384000 / 2)) {
						param.nr = temp;
					}
				}
				if (sscanf(p1,"-nr_option:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 5) {
						param.nrLV = (int)temp - 1;
					}
				}
				if (sscanf(p1,"-hfa_preset:%d",&temp) == 1) {
					if (temp >= 1 && temp <= 10) {
						param.hfa_preset = temp;
					}
				}

				if (sscanf(p1,"-thread:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 24) {
						param.thread = (int)temp;
					}
				}

				if (strcmp(p1,"-hfaFast") == 0) {
					param.hfaFast = 1;
				}
				if (strcmp(p1,"-hfaWide") == 0) {
					param.hfaWide = 1;
				}
				if (sscanf(p1,"-deE:%ld",&temp) == 1) {
					if (temp == 0 || temp == 1 || temp == 2) {
						param.deEmphasis = (int)temp;
					}
				}
				if (strcmp(p1,"-C") == 0 || strcmp(p1,"-SLR") == 0 || strcmp(p1,"-LFE") == 0) {
					param.mc_flag = 1;
				}

#if 0
				if (strcmp(p1,"-adjFreq") == 0) {
					param.ana = 1;
				}
				if (strcmp(p1,"-spAna") == 0) {
					// 192khz 固定で解析
					param.sp_ana = 1;
				}
#endif
				if (sscanf(p1,"-fio:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 16000) {
						param.fio = temp;
					}
				}
				if (sscanf(p1,"-cpu_pri:%d",&temp) == 1 && temp >= 0 && temp <= 1) {
					param.cpu_pri = temp;
				}
//				if (strcmp(p1,"-eq") == 0) {
//					param.eq_flag = 1;
//				}
//				if (sscanf(p1,"-temp:%c%c%c",tmppath,tmppath+1,tmppath+2) == 3) {
//					param.tempPath[0] = tmppath[0];
//					param.tempPath[1] = tmppath[1];
//					param.tempPath[2] = '\0';
//					param.tempPath[3] = '\0';
//				}
//				if (sscanf(p1,"-dsd_fmt:%ld",&temp) == 1) {
//					if ((temp == 64 || temp == 128 || temp == 256)) {
//						param.dsd_fmt = temp;
//					}
//				}

				if (sscanf(p1,"-upconv:%ld",&temp) == 1) {
					param.upconv = temp;
				}
				
				if (p2 == NULL) {
					break;
				}
				p1 = p2 + 1;
				p2 = strchr(p1,(int)' ');
			}
			if (param.enable_hfc == 0) {
				param.hfc = -1;
			}
			if (param.enable_lfc == 0) {
				param.lfc = -1;
			}
			if (param.enable_nr == 0) {
				param.nr = -1;
			}
#if 0
			if (param.sp_ana == 1) {
				param.outSampleR = 192000;
				param.hfa		 = 0;
				param.overSamp	 = 0;
				param.hfc		 = -1;
				param.lfc		 = -1;
				param.lpf		 = -1;
				param.nr		 = -1;
				param.hfc_auto	 = 0;
				param.eq_flag	 = 0;
				param.post_abe	 = 0;
			}
#endif
			if (param.upconv == 0) {
				int ip;
				int retStatus;
				char *arg[5];
				char buffer[2048];
				// wav2raw 変換処理
				paramFlag = 0x07;
				arg[0] = strdup(argv[0]);
				arg[1] = malloc(4096);if (arg[1] != NULL) strcpy(arg[1],argv[1]);
				arg[2] = strdup(argv[2]);
				arg[3] = strdup(argv[3]);
				arg[4] = (char *)malloc(2048);
				if (arg[0] == NULL || arg[1] == NULL || arg[2] == NULL || arg[3] == NULL || arg[4] == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param.errLine = __LINE__;
					break;
				}
				strcpy(arg[4],argv[4]);

				retStatus = to_raw_main(5,arg);
				if (retStatus != STATUS_SUCCESS) {
					break;
				}
				p1 = strstr(arg[4],"-ch:");
				if (p1 != NULL) {
					if (sscanf(p1,"-ch:%d",&temp) == 1) {
						param.channel_count = temp;
					}
				}
				if (param.channel_count == 0) {
					retCode = STATUS_UNKNOWN_FORMAT;param.errLine = __LINE__;
					break;
				}

				for (ip = 0;ip < param.channel_count;ip++) {
					// 子プロセス
					char param4[4096];
					char s[20];
					sprintf(s," -upconv:%d",ip + 1);
					p1 = strstr(arg[4],"-upconv:");
					if (p1 != NULL) {
						sprintf(p1,"-upconv:%d",ip + 1);
					} else {
						strcat(arg[4],s);
					}
					if (param.mc_flag) {
						sprintf(s," -ms:%d",param.outSampleR * 2);
						strcat(arg[4],s);
					}

					if (start_exec(5,arg,param.cpu_pri,&startexec_info[ip].hStdOutRead,&startexec_info[ip].ps,&startexec_info[ip].thread) == -1) {
						retCode = STATUS_EXEC_FAIL;param.errLine = __LINE__;
						break;
					}
					startexec_info[ip].state = 1;
				}
				if (retCode != 0) break;

				for (;;) {
					unsigned long read_size;
					if (!ReadFile(startexec_info[0].hStdOutRead,buffer,sizeof(buffer),&read_size,NULL)) {
						if (GetLastError() == ERROR_BROKEN_PIPE) {
							break; // pipe done
						} else {
							// Error
							break;
						}
					}
					if (read_size == 0) {
						continue;
					}
					buffer[read_size] = '\0';
					fprintf(stdout,"%s",buffer);
					fflush(stdout);
					for (ip = 1;ip < param.channel_count;ip++) {
						ReadFile(startexec_info[ip].hStdOutRead,buffer,sizeof(buffer),&read_size,NULL);
					}
				}
				if (1) {
					HANDLE ph[6];
					int all;
					all = 1;
					for (ip = 0;ip < param.channel_count;ip++) {
						if (startexec_info[ip].state == 1) {
							ph[ip] = startexec_info[ip].ps;
						} else {
							all = 0;
							break;
						}
					}
					
					if (all) {
						WaitForMultipleObjects(param.channel_count,ph,1,0xFFFFFFFF);
						for (ip = 0;ip < param.channel_count;ip++) {
							finish_exec(&startexec_info[ip].hStdOutRead,&startexec_info[ip].ps,&startexec_info[ip].thread);
						}
					} else {
						for (ip = 0;ip < param.channel_count;ip++) {
							if (startexec_info[ip].state == 1) {
								WaitForSingleObject(startexec_info[ip].ps,0xFFFFFFFF);
								finish_exec(&startexec_info[ip].hStdOutRead,&startexec_info[ip].ps,&startexec_info[ip].thread);
							}
						}
						retCode = STATUS_EXEC_FAIL;param.errLine = __LINE__;
						break;
					}
				}
				param.upconv = 0;
				if (param.mc_flag == 1) {
					retStatus = mc_main(5,arg);
					if (retStatus != STATUS_SUCCESS) {
						break;
					}
				}
				retStatus = to_wav_main(5,arg);
				if (retStatus != STATUS_SUCCESS) {
					break;
				}
				
			} else {
				// r1 ファイルかをチェックする(sp_ana)
				if (param.upconv == 1) {
					param.r1_flag = 1;
				}

				// GenNoise 用のデータパス
				_splitpath(argv[0],drive,dir,fname,ext);
				_makepath(tmppath,drive,dir,"nd","dat");
				param.nd_path = malloc(strlen(tmppath) + 1);
				if (param.nd_path != NULL) {
					strcpy(param.nd_path,tmppath);
				}

				// ビット拡張テーブル
				_splitpath(argv[0],drive,dir,fname,ext);
				_makepath(tmppath,drive,dir,"bit_extend_table","");
				
				fp = fopen(tmppath,"r");
				if (fp) {
					for (i = 0;i < 256;i++) {
						if (fscanf(fp,"%d,",&c) == 1) {
							param.beInfo.cutOff[i] = c;
						}
					}
					fclose(fp);
				}

				// hfa2/3 Preset
				_splitpath(argv[3],drive,dir,fname,ext);
				_makepath(tmppath,drive,dir,"preset","dat");
				fp = fopen(tmppath,"r");
				if (fp) {
					for (i = 1;i <= param.hfa_preset;i++) {
						if (fgets(work,2000,fp) != NULL) {
							if (sscanf(work,"%[^,],%d,%d,%d,%d,%d,%d,%d,%d,%d",workstr,&temp,&temp2,&temp3,&temp4,&temp5,&temp6,&temp7,&temp8,&temp9) == 10) {
								if (temp == 0) {
									param.sig1Enb = 0;
								} else {
									param.sig1Enb = 1;
								}
								if (temp2 >= 1 && temp2 <= 25) param.sig1AvgLineNx = temp2;
								if (temp3 >= -44 && temp3 <= 44) param.sig1Phase = temp3;
								if (temp4 == 0) {
									param.sig2Enb = 0;
								} else {
									param.sig2Enb = 1;
								}
								if (temp5 >= -44 && temp5 <= 44) param.sig2Phase = temp5;
								if (temp7 >= 0 && temp7 <= 100)  param.hfaNB = temp7;
								if (temp8 == 0) {
									param.sig3Enb = 0;
								} else {
									param.sig3Enb = 1;
								}
								if (temp9 >= -44 && temp9 <= 44) param.sig3Phase = temp9;
							}
						}
					}
					fclose(fp);
				}
				sprintf(work,"HFA Preset:%d,%d,%d,%d,%d,%d,%d,%d",param.sig1Enb,param.sig1AvgLineNx,param.sig1Phase,param.sig2Enb,param.sig2Phase,param.sig3Enb,param.sig3Phase,param.hfaNB);
				PRINT_LOG(work);
#ifdef _OPENMP
				int nCpu;
				nCpu = param.thread;
				omp_set_num_threads(nCpu);
#endif
				if (param.hfc != -1) {
					param.hfc_auto = 0;
				}

				_splitpath(argv[2],drive,dir,fname,ext);
				_makepath(tmppath,drive,dir,fname,"abort");
				strcpy(param.abort_filename,tmppath);

				_splitpath(workpath,workdrive,workdir,workfname,workext);
				_splitpath(argv[2],drive,dir,fname,ext);
				if (strlen(workpath) >= 3) {
					_makepath(workpath,workdrive,workdir,fname,ext);
					rn_fname = workpath;
				} else {
					rn_fname = argv[2];
				}

				// 入力ファイル名の生成
				_splitpath(rn_fname,drive,dir,fname,ext);
				sprintf(ext,"r%d",param.upconv);
				_makepath(rn_fname,drive,dir,fname,ext);

				if (paramFlag == 0x0007) {
					SamplingRateConvert(rn_fname,&param);
					if (param.err == STATUS_SUCCESS) {
						if (NormInfo.max < 0) {
							NormInfo.max *= -1;
						}
						if (NormInfo.min < 0) {
							NormInfo.min *= -1;
						}
						max = NormInfo.max;
						if (max < NormInfo.min) {
							max = NormInfo.min;
						}
						percent = (double)max / (double)0x7FFFFFFFFFFFFF;
						avg = NormInfo.avg << 40;
						strcpy(tmppath,rn_fname);
						strcat(tmppath,".param");
						fp = fopen(tmppath,"wb");
						if (fp == NULL) {
							param.err = STATUS_FILE_WRITE_ERR;
							break;
						}
						if (fprintf(fp,"%.10lf,%llx\n",percent,avg) == EOF) {
							param.err = STATUS_FILE_WRITE_ERR;
							break;
						}
						fclose(fp);
					}
				}
			}
		}
		if (argc != 5 || paramFlag != 0x0007) {
			printf(STR_COPYRIGHT);
			printf(STR_USAGE);
			exit(0);
		}
	} while (0);

	if (param.upconv != 0 && param.err != STATUS_SUCCESS) {
		_splitpath(argv[2],drive,dir,fname,ext);
		_makepath(tmppath,drive,dir,fname,"err");
		fp = fopen(tmppath,"a");
		if (fp) {
			switch (param.err) {
				case STATUS_PARAMETER_ERR:
					fprintf(fp,"%s[%d]:Parameter error.\n",param.upconv == 1 ? "upconv" : "wav2wav",param.errLine);
					break;
				case STATUS_FILE_READ_ERR:
					fprintf(fp,"%s[%d]:File read error.\n",param.upconv == 1 ? "upconv" : "wav2wav",param.errLine);
					break;
				case STATUS_FILE_WRITE_ERR:
					fprintf(fp,"%s[%d]:File write error.\n",param.upconv == 1 ? "upconv" : "wav2wav",param.errLine);
					break;
				case STATUS_MEM_ALLOC_ERR:
					fprintf(fp,"%s[%d]:Memory Allocation error.\n",param.upconv == 1 ? "upconv" : "wav2wav",param.errLine);
					break;
				default:
					fprintf(fp,"%s[%d]:Other error.\n",param.upconv == 1 ? "upconv" : "wav2wav",param.errLine);
			}
			fclose(fp);
		}
	}
	if (param.upconv == 0) {
		_splitpath(argv[2],drive,dir,fname,ext);
		_makepath(tmppath,drive,dir,fname,"files");
		fp_files = fopen(tmppath,"r");
		if (fp_files != NULL) {
			while (fgets(workpath,_MAX_PATH,fp_files) != NULL) {
				p1 = strrchr(workpath,'\n');if (p1 != NULL) *p1 = '\0';
//				fprintf(stdout,"Delete:%s\n",workpath);
//				fflush(stdout);
				PRINT_LOG(workpath);
				unlink(workpath);
			}
			fclose(fp_files);
			_splitpath(argv[2],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,fname,"files");
			PRINT_LOG(tmppath);
			unlink(tmppath);
		}
		fprintf(stdout,"\n");
		fflush(stdout);
		fprintf(stdout,"End\n");
		fflush(stdout);
	}
	exit(retCode);

	return 0;
}
//---------------------------------------------------------------------------
int is_ext_audiofile(char *filename,char *ext)
{
	char workdrive[_MAX_DRIVE];
	char workdir[_MAX_DIR];
	char workfname[_MAX_FNAME];
	char workext[_MAX_EXT];
	
}

//---------------------------------------------------------------------------
// Function   : SamplingRateConvert
// Description: サンプリングレート変換処理をする
// ---
//	rawFile	: RAWデータファイル名
//	param	: 変換パラメータ構造体
//
void SamplingRateConvert(char *rawFile,PARAM_INFO *param)
/*
 *	サンプリングレート変換
 */
{
	char outFile[_MAX_PATH];
	char tmpFile[_MAX_PATH];
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];
	FIO fp_r;
	FIO fp_w;
	FIO fp_r2;
	fio_size size;
	SSIZE inSample,outSample;
	int hfc;
	long double dd;
	DWORD svInSampleR;
	int svIw;
	int downSampleFlag = FALSE;
	DWORD svOutSampleR;
	int validIndex;
	int i;
	FFT_PARAM fftParam;

	memset(&fp_r,0,sizeof (FIO));
	memset(&fp_r2,0,sizeof (FIO));
	memset(&fp_w,0,sizeof (FIO));

	memset(&fftParam,0,sizeof (FFT_PARAM));
	fftParam.inSampleR  = param->inSampleR;
	fftParam.outSampleR = param->outSampleR;
	fftParam.hfc        = param->hfc;
	fftParam.lfc        = param->lfc;
	fftParam.cut_high_dither   = param->cut_high_dither;

	fftParam.eq_ref_max = (double *)malloc(192001 * sizeof (double));
	fftParam.eq_ref_avg = (double *)malloc(192001 * sizeof (double));
	fftParam.eq         = (double *)malloc(192001 * sizeof (double));
	fftParam.lfa_eq     = (double *)malloc(192001 * sizeof (double));
	fftParam.eq_ref_count = 0;
	fftParam.deEmphasis = param->deEmphasis;
	fftParam.dsd_fmt    = -1;
	fftParam.abort_filename = param->abort_filename;
	fftParam.abort_percent  = 0;

	do {
		if (fftParam.eq_ref_max == NULL || fftParam.eq_ref_avg == NULL || fftParam.eq == NULL || fftParam.lfa_eq == NULL) {
			param->err = STATUS_MEM_ALLOC_ERR;
			break;
		}

		// テンポラリファイル名の作成
//		if (strlen(param->tempPath) > 0) {
//			_splitpath(rawFile,drive,dir,fname,ext);
//			strcpy(dir,"\\");
//			_makepath(outFile,param->tempPath,dir,fname,ext);
//		} else {
			strcpy(outFile,rawFile);
//		}
		strcat(outFile,".tmp");
		strcpy(tmpFile,outFile);
		strcat(tmpFile,"2");
		
		// リード専用ファイルオープン関数(バッファリング機能付)
		fio_open(&fp_r,rawFile,FIO_MODE_R);
		if (fp_r.error) {
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		if (param->dsd_fmt == -1) {
			if (param->overSamp == 1) {
				// 2倍のサンプリングレートで計算し後でダウンサンプルする。
				svOutSampleR = fftParam.outSampleR;
				fftParam.outSampleR *= 2;
				param->outSampleR *= 2;
			} else if (param->overSamp == 2) {
				// 768000*2のサンプリングレートで計算し後でダウンサンプルする。
				svOutSampleR = fftParam.outSampleR;
				fftParam.outSampleR = 384000;
				param->outSampleR   = 384000;
			} else if (param->overSamp == 3) {
				svOutSampleR = fftParam.outSampleR;
				fftParam.outSampleR = 1536000;
				param->outSampleR   = 1536000;
			} else if (param->overSamp == 4) {
				// 768000*2のサンプリングレートで計算し後でダウンサンプルする。
				svOutSampleR = fftParam.outSampleR;
				fftParam.outSampleR = 2822400;
				param->outSampleR   = 2822400;
			}
		}
		outSample = 0;
		
		// ファイルサイズ取得
		fio_get_filesize(&fp_r,&size);
		if (fp_r.error) {
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		
		// 出力サンプル数を計算する。
		inSample = (SSIZE)(size / sizeof (SSIZE));
		dd = inSample;
		dd *= fftParam.outSampleR;
		dd /= fftParam.inSampleR;
		outSample = (SSIZE)dd;

		if (outSample == 0) {
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		// 量子化ビット拡張処理
		if (param->abe != 0) {
			fprintf(stdout,"[ABE]\n");
			fflush(stdout);
			// 出力用にファイルオープン
			fio_open(&fp_w,outFile,FIO_MODE_W);
			if (fp_w.error) {
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			fio_set_memory_limit(&fp_w,20,param->fio);
			adjBitExtension(inSample,&fp_r,&fp_w,param);
			if (param->err) {
				break;
			}
			fio_close(&fp_r);
			if (fp_r.error) {
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
			fio_setmode_r(&fp_w,&fp_r,rawFile);
			if (fp_w.error) {
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			if (fp_r.error) {
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
		}

		fprintf(stdout,"[SRC]\n");
		fflush(stdout);
		if (inSample > outSample) {
			downSampleFlag = TRUE;
		}

		fio_open(&fp_w,outFile,FIO_MODE_W);
		if (fp_w.error) {
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}

		// ファイルに出力するサイズを制限する(outSample数)
		fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
		fio_set_memory_limit(&fp_w,20,param->fio);

		for (i = 0;i < 192000;i++) {
			fftParam.eq_ref_max[i] = 0;
			fftParam.eq_ref_avg[i] = 0;
			fftParam.eq[i] = 0;
			fftParam.lfa_eq[i] = 0;
		}

		fftParam.analyze_mode = 1;
		if (param->overSamp >= 3) {
			fftParam.hi_sample_mode = 1;
		}
		fftFilter(inSample,outSample,&fp_r,&fp_w,&fftParam);
		if (param->err) {
			break;
		}

		if (fftParam.analyze_mode == 1) {
			int i,j,n;
			double avg1,avg2,p1_avg,p2_avg,p1_step,p1_per,p2_per;
			int p1_index,p2_index,p1_e_index,p2_e_index;

			validIndex = param->inSampleR / 2;
			if (validIndex > 192000) validIndex = 192000;

			if (fftParam.eq_ref_count > 0) {
				for (i = 1;i < 192000;i++) {
					fftParam.eq_ref_avg[i] /= fftParam.eq_ref_count;
				}
			}
			for (i = 0;i < validIndex - 10;i++) {
				avg1 = 0;
				avg2 = 0;
				for (j = 0;j < 10;j++) {
					avg1 += fftParam.eq_ref_max[i + j];
					avg2 += fftParam.eq_ref_avg[i + j];
				}
				avg1 /= 10;
				avg2 /= 10;
				for (j = 0;j < 10;j++) {
					fftParam.eq_ref_max[i + j] = avg1;
					fftParam.eq_ref_avg[i + j] = avg2;
				}
			}

			if (0) {
				FILE *ff;
				ff = fopen("d:\\ref_max.txt","wt");
				if (ff != NULL) {
					for (i = 0;i < 192000;i++) {
						fprintf(ff,"%lf\n",fftParam.eq_ref_max[i]);
					}
					fclose(ff);
				}
				ff = fopen("d:\\ref_avg.txt","wt");
				if (ff != NULL) {
					for (i = 0;i < 192000;i++) {
						fprintf(ff,"%lf\n",fftParam.eq_ref_avg[i]);
					}
					fclose(ff);
				}
			}

			if (fftParam.cut_high_dither == 1) {
char sssss[256];
				p1_index   = 5000;
				p1_e_index = 6000;
				p2_index   = 10000;
				p2_e_index = 11000;
				p1_avg = 0;
				p2_avg = 0;
				for (i = p1_index,n = 0;i < p1_e_index;i++,n++) {
					p1_avg += fftParam.eq_ref_avg[i];
				}
				if (n > 0) {
					p1_avg /= n;
				}
				for (i = p2_index,n = 0;i < p2_e_index;i++,n++) {
					p2_avg += fftParam.eq_ref_avg[i];
				}
				if (n > 0) {
					p2_avg /= n;
				}

				// 4kHz ～ 5kHz の間のpowerと8kHz～9kHzのpowerを比べ、高音が大きい場合は調整用のデータを生成する。
				p1_per = p2_avg / p1_avg;
				if (p1_per > 0.53) {
					p1_index   = 6000;
					p2_index   = 10000;
					p1_per = ((double)0.47 / p1_per);
					p1_step = ((double)1.0 - p1_per) / (p2_index - p1_index);
					for (i = 1;i < validIndex;i++) {
						if (i > p2_index) {
							fftParam.eq[i] = fftParam.eq[i - 1];
						} else if (i > p1_index) {
							fftParam.eq[i] = fftParam.eq[i - 1] - p1_step;
						} else {
							fftParam.eq[i] = 1;
						}
					}
				} else {
					fftParam.cut_high_dither = 0;
				}
			}
		}

		// 音量調査
		fio_flush(&fp_w);
		merageTempFile(' ',1,&fp_w,NULL,NULL,outSample,param);
		if (param->err) {
			break;
		}

		fio_close(&fp_w);

		// hfc auto 用パラメーター設定(hfc)
		if (param->hfc_auto == 1) {
			anaHFC_AutoParam(&fftParam);
			param->hfc = fftParam.hfc;
		}
		// lfa 用パラメーター作成(lfa_eq)
		if (param->low_adjust != 0) {
			anaLFA_Param(&fftParam);
		}

		fio_open(&fp_w,outFile,FIO_MODE_W);
		if (fp_w.error) {
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}
		// ファイルに出力するサイズを制限する(outSample数)
		fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
		fio_set_memory_limit(&fp_w,20,param->fio);

		// fft
		fftParam.analyze_mode = 0;
		fftParam.lfa_flag = param->low_adjust;
		fftParam.lpf_flag = param->lpf;
		fftFilter(inSample,outSample,&fp_r,&fp_w,&fftParam);
		if (param->err) {
			break;
		}
		fio_close(&fp_r);
		if (fp_r.error) {
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		

		// 生成したデータを読み込み用に設定する。
		fio_setmode_r(&fp_w,&fp_r,rawFile);
		if (fp_w.error) {
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}
		if (fp_r.error) {
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		param->cut_high_dither = 0;
		param->low_adjust = 0;
		param->lpf = 0;
		param->ana = 0;
		param->eq_flag = 0;
		fftParam.lfa_flag = 0;
		fftParam.lpf_flag = 0;
		fftParam.cut_high_dither = 0;

//		if (param->eq_pw != NULL) {
//			free(param->eq_pw);
//		}

		// スピーカー用周波数解析
//		if (param->sp_ana == 1) {
//			spAnalyze(outSample,&fp_r,param);
//		}

		//
		// ノイズカット
		if (param->nr != -1) {
			fprintf(stdout,"[NR]\n");
			fflush(stdout);
			if (param->hfc != -1) {
				hfc = param->hfc;
			} else {
				if (downSampleFlag == FALSE) {
					hfc = param->inSampleR / 2;
				} else {
					hfc = param->outSampleR / 2;
				}
			}
			if (downSampleFlag == FALSE) {
				if (hfc > param->inSampleR / 2) {
					hfc = param->inSampleR / 2;
				}
			} else {
				if (hfc > param->outSampleR / 2) {
					hfc = param->outSampleR / 2;
				}
			}
			if (param->nr < hfc) {
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,20,param->fio);

				noiseCut(param->nr,outSample,&fp_r,&fp_w,param);
				if (param->err) {
					break;
				}
				fio_setmode_r(&fp_w,&fp_r2,tmpFile);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r2.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,20,param->fio);

				//bpFilter(-1,hfc,outSample,&fp_r2,&fp_w,param);
				fftParam.disable_eq = 1;
				fftParam.lfc = -1;
				fftParam.hfc = hfc;
				fftParam.inSampleR = fftParam.outSampleR;
				fftFilter(outSample,outSample,&fp_r2,&fp_w,&fftParam);
				if (param->err) {
					break;
				}
				fio_close(&fp_r2);
				fio_setmode_r(&fp_w,&fp_r2,tmpFile);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r2.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,20,param->fio);

				merageTempFile('-',0,&fp_r,&fp_r2,&fp_w,outSample,param);
				if (param->err) {
					break;
				}
				fio_close(&fp_r);
				fio_close(&fp_r2);
				if (fp_r.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				if (fp_r2.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_setmode_r(&fp_w,&fp_r,rawFile);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
			}
		}

		//
		// 高域補間処理
		if (param->hfa != 0) {
			fflush(stdout);

			if (downSampleFlag == FALSE) {
				hfc = param->inSampleR / 2;
				if (param->hfc != -1 && hfc > param->hfc) {
					hfc = param->hfc;
				}
			} else {
				hfc = param->outSampleR / 2;
				if (param->hfc != -1 && hfc > param->hfc) {
					hfc = param->hfc;
				}
			}

			if (param->hfa != 1) {
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,20,param->fio);

				fprintf(stdout,"[HFA%d]\n",param->hfa);
				fflush(stdout);
				fftParam.disable_eq = 1;
				fftParam.lfc = -1;
				fftParam.hfc = hfc;
				fftParam.inSampleR = fftParam.outSampleR;
				genOverTone(hfc,outSample,&fp_r,&fp_w,param,&fftParam);
				if (param->err) {
					break;
				}
				fio_setmode_r(&fp_w,&fp_r2,tmpFile);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r2.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,20,param->fio);

				//bpFilter(hfc,-1,outSample,&fp_r2,&fp_w,param);
				fftParam.disable_eq = 1;
				fftParam.lfc = hfc;
				fftParam.hfc = -1;
				fftParam.inSampleR = fftParam.outSampleR;
				fftFilter(outSample,outSample,&fp_r2,&fp_w,&fftParam);
				if (param->err) {
					break;
				}
				fio_close(&fp_r2);
				fio_setmode_r(&fp_w,&fp_r2,tmpFile);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r2.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,20,param->fio);

				merageTempFile('+',0,&fp_r,&fp_r2,&fp_w,outSample,param);
				if (param->err) {
					break;
				}
				fio_close(&fp_r);
				fio_close(&fp_r2);
				if (fp_r.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				if (fp_r2.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_setmode_r(&fp_w,&fp_r,rawFile);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
			}
			if (param->hfa == 1 || param->hfaNB > 0) {
				fprintf(stdout,"[HFA1]\n");
				fflush(stdout);
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,20,param->fio);

				genNoise(hfc,outSample,&fp_r,&fp_w,param);
				if (param->err) {
					break;
				}
				fio_close(&fp_r);
				if (fp_r.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_setmode_r(&fp_w,&fp_r,rawFile);
				if (fp_w.error) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r.error) {
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
			}
		}
		if (param->post_abe == 1) {
			fprintf(stdout,"[Post ABE]\n");
			fflush(stdout);
			fio_open(&fp_w,outFile,FIO_MODE_W);
			if (fp_w.error) {
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			// ファイルに出力するサイズを制限する(outSample数)
			fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
			fio_set_memory_limit(&fp_w,20,param->fio);

			svInSampleR = param->inSampleR;
			svIw		= param->iw;
			param->inSampleR = param->outSampleR;
			param->iw = param->w;
			param->cutLowData = 0;
			param->abeFnLevel = 0;
			param->smLowData  = 0;
			param->adaptiveFilder = 0;
			param->abeNX = 0;
			adjBitExtension(outSample,&fp_r,&fp_w,param);
			if (param->err) {
				break;
			}
			param->inSampleR = svInSampleR;
			param->iw = svIw;

			fio_close(&fp_r);
			fio_setmode_r(&fp_w,&fp_r,rawFile);
			if (fp_w.error) {
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			if (fp_r.error) {
				param->err = STATUS_FILE_READ_ERR;
				break;
			}

		}

		// オーバーサンプリング用のダウンサンプル処理
		if (param->dsd_fmt == -1 && param->overSamp != 0) {
			char sssss[128];
			SSIZE wkSample;
			memset(&NormInfo,0,sizeof (NORM_INFO));

			fio_open(&fp_w,outFile,FIO_MODE_W);
			if (fp_w.error) {
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			fio_set_memory_limit(&fp_w,20,param->fio);

			fprintf(stdout,"[SRC(DS)]\n");
			fflush(stdout);

			fftParam.hfc = -1;
			fftParam.lfc = -1;
			fftParam.hi_sample_mode = 1;
			fftParam.inSampleR  = param->outSampleR;
			param->inSampleR    = param->outSampleR;
			fftParam.outSampleR = svOutSampleR;
			param->outSampleR   = svOutSampleR;
			dd = outSample;
			dd *= param->outSampleR;
			dd /= param->inSampleR;
			wkSample = (SSIZE)dd;
			// ファイルに出力するサイズを制限する(outSample数)
			fio_set_maxsize(&fp_w,(fio_size)wkSample * sizeof (SSIZE));

			fftParam.disable_eq = 1;
			//fftParam.lvadj_flag = 1;
			fftFilter(outSample,wkSample,&fp_r,&fp_w,&fftParam);
			if (param->err) {
				break;
			}
			fio_close(&fp_r);
			if (fp_r.error) {
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
			fio_setmode_r(&fp_w,&fp_r,rawFile);
			if (fp_w.error) {
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			if (fp_r.error) {
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
			merageTempFile(' ',1,&fp_r,NULL,NULL,wkSample,param);
			if (param->err) {
				break;
			}
		} else if (param->dsd_fmt != -1) {
			SSIZE wkSample;
			memset(&NormInfo,0,sizeof (NORM_INFO));

			fio_open(&fp_w,outFile,FIO_MODE_W);
			if (fp_w.error) {
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			fio_set_memory_limit(&fp_w,20,param->fio);

			fprintf(stdout,"SRC(DSD)\n");
			fflush(stdout);

			param->inSampleR    = param->outSampleR;
			fftParam.inSampleR  = param->outSampleR;
			if (param->dsd_fmt == 64) {
				fftParam.hfc = 34100;
				param->outSampleR   = 2822400;
				fftParam.outSampleR = 2822400;
			} else if (param->dsd_fmt == 128) {
				fftParam.hfc = 44100;
				param->outSampleR   = 2822400 * 2;
				fftParam.outSampleR = 2822400 * 2;
			} else if (param->dsd_fmt == 256) {
				fftParam.hfc = 66100;
				param->outSampleR   = 2822400 * 4;
				fftParam.outSampleR = 2822400 * 4;
			}
			fftParam.lfc = 20;
			fftParam.hi_sample_mode = 1;
			fftParam.dsd_fmt = param->dsd_fmt;
			fftParam.analyze_mode = 0;
			dd = outSample;
			dd *= param->outSampleR;
			dd /= param->inSampleR;
			wkSample = (SSIZE)dd;
			// ファイルに出力するサイズを制限する(outSample数)
			fio_set_maxsize(&fp_w,(fio_size)wkSample * sizeof (SSIZE));

			fftParam.disable_eq = 1;
			fftFilter(outSample,wkSample,&fp_r,&fp_w,&fftParam);
			if (param->err) {
				break;
			}
			fio_close(&fp_r);
			if (fp_r.error) {
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
			fio_setmode_r(&fp_w,&fp_r,rawFile);
			if (fp_w.error) {
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			if (fp_r.error) {
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
			merageTempFile(' ',1,&fp_r,NULL,NULL,wkSample,param);
			if (param->err) {
				break;
			}
		}
		fio_close(&fp_r);
		//
		// 終わり
	} while (0);
}

//---------------------------------------------------------------------------
// Function   : anaLfaParam
// Description: LFA パラメーター作成
// ---
//	param		:変換パラメータ
//
void anaLFA_Param(FFT_PARAM *param)
{
	int i,ii,n;
	long double ref_pw,pw;

	ref_pw = 0;
	for (i = 200,n = 0;n < 20;i++,n++) {
		ref_pw += param->eq_ref_avg[i];
	}
	if (n > 0) {
		ref_pw /= n;
	}

	for (i = 210;i >= 20;) {
		pw = 0;
		for (ii = i,n = 0;n < 20;ii++,n++) {
			pw += param->eq_ref_avg[ii];
		}
		if (n > 0) {
			pw /= n;
		}
		if ((pw / ref_pw) >= 0.8 && (pw / ref_pw) <= 1.0) {
			for (ii = i,n = 0;n < 20;ii++,n++) {
				param->lfa_eq[ii] = 0.96 / (pw / ref_pw);
			}
		}
		i -= n;
	}
}
//---------------------------------------------------------------------------
// Function   : anaHFC_autoParam
// Description: HFC Auto パラメーター作成
// ---
//	param		:変換パラメータ
//
void anaHFC_AutoParam(FFT_PARAM *param)
{
	
	long index,i,n;
	double p;
	static double pw[192001];

	for (i = 0;i < 192000;i++) {
		pw[i]  = 0;
	}

	// 200Hz ごとの平均をとる
	for (index = 1;index + 100 < 192000;index += 100) {
		p = 0;
		for (i = index,n = 0;n < 100;n++,i++) {
			p += param->eq_ref_avg[i];
		}
		if (n > 0) {
			p /= n;
		}
		for (i = index,n = 0;n < 100;n++,i++) {
			pw[i] = p;
		}
	}

	for (i = 192000;i > 0;i--) {
		if (pw[i] < (double)10000000000000000) {
			pw[i] = 0;
		} else {
			break;
		}
	}

	for (i = 192000;i > 0;i--) {
		if (pw[i] != 0) {
			break;
		}
	}

	if (i < param->inSampleR / 2) {
		param->hfc = i;
		if (param->hfc < 11000) param->hfc = 11000;
	} else {
		param->hfc = param->inSampleR / 2;
	}
}

//---------------------------------------------------------------------------
// Function   : adjBitExtension
// Description: ビット分解能を高める処理
// ---
//	inSample 	:処理するサンプル数
//	fp			:入力ファイル
//	tmpFp		:出力ファイル
//	param		:変換パラメータ
//
void adjBitExtension(SSIZE inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3,*mem4;
	SSIZE level,level2;
	double nx;
	long wkMemSize;
	long inSampleR;
	long fftSize,i,j,k,n;
	long cnt;
	long cutOff;
	long hfc;
	long wr;
	double percent,per;
	SSIZE *pIn,*pIn2,*pOut;
	SSIZE d1,d2,d3,d4,d5,d6,d7,d8;
	SSIZE dd1,dd2,dd3,dd4;
	SSIZE startInSample,nSample;
	SSIZE outRemain;
	double samp[9];
	double power;
	double shresh;
	int next;
	fftw_complex *fftw_in,*fftw_out;
	fftw_plan fftw_p,fftw_ip;
	int sm_ignore;
	SSIZE sm_avg;

	fio_rewind(fp_r);

	inSampleR = param->inSampleR;

	fftSize = inSampleR * 2;
	wkMemSize = fftSize * 2;

	mem1 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem2 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem3 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem3 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem4 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem4 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_in == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_out == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p = fftw_plan_dft_1d(fftSize,fftw_in,fftw_out,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip = fftw_plan_dft_1d(fftSize,fftw_out,fftw_in,FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	if (param->iw == 16) {
		shresh = (double)12592500000000000.0;	// 1k
		shresh = (double)4200000000000000.0;	// 1k
	} else if (param->iw == 24) {
		shresh = (double)85120000000000.0;	// 1k
	} else {
		shresh = -1;
	}
	outRemain = inSample;
	per = -1;
	for (startInSample = (((fftSize * 2) + (fftSize / 2)) * -1);startInSample < inSample + (fftSize * 1);startInSample += fftSize) {
		if (startInSample >= 0 && startInSample  < inSample) {
			percent = ((double)startInSample / inSample);
			percent *= 100;
			percent = (double)((int)percent);
			if (percent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)percent);
				fflush(stdout);
				if (chkAbort(param,percent,10) == 1) exit(0);
			}
			per = percent;
		}

		memset(mem1,00,wkMemSize * sizeof (SSIZE));
		nSample = fftSize * 2;

		if (startInSample >= 0 && startInSample + (fftSize * 1) < inSample + fftSize) {
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			nSample = fftSize * 2;
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else if (startInSample + (fftSize * 1) >= 0 || startInSample < inSample) {
			mem0 = mem1;
			nSample = 0;
			if (startInSample < 0) {
				fio_seek(fp_r,0,SEEK_SET);
				if ((startInSample * -1) < fftSize * 2) {
					mem0 += (startInSample * -1);
					nSample = (fftSize * 2) + startInSample;
				}
			} else if (startInSample < inSample) {
				fio_seek(fp_r,startInSample,SEEK_SET);
			}
			if (nSample > 0 && startInSample < inSample && startInSample + nSample > inSample) {
				nSample = inSample - startInSample;
			}

			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}
		
		// 音のレベルを調査しておく
		if (param->abeNX) {
			level = 0;
			for (i = fftSize / 2,j = 0,n = 0;n < fftSize;i++,j++,n++) {
				if (mem1[i] > 0) {
					level += mem1[i] >> (56 - 16);
					j++;
				}
			}
			if (j > 0) {
				level /= j;
			}
		}
#if 1
		if (param->abeFnLevel > 0) {
			// ディザキャンセラ
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			// ２値平均
			for (i = 0;i + 2 < nSample;i++) {
				d1 = pIn[i + 0];
				d2 = pIn[i + 1];
				d3 = pIn[i + 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
//				if (dd1 < 0) {
//					dd1 *= -1;
//				}
//				if (dd2 < 0) {
//					dd2 *= -1;
//				}
//				if (dd3 < 0) {
//					dd3 *= -1;
//				}
				dd1 /= param->abeFnLevel;
				dd2 /= param->abeFnLevel;
				dd3 /= param->abeFnLevel;
				sm_ignore = 1;
				if (dd1 + 1 == dd2 && dd2 == dd3 + 1) {
					sm_ignore = 0;
				}
				if (dd1 == dd2 + 1 && dd2 + 1 == dd3) {
					sm_ignore = 0;
				}

				if (sm_ignore == 0) {
					sm_avg = 0;
					sm_avg += pIn[i + 0];
					sm_avg += pIn[i + 1];
					sm_avg /= 2;
					pOut[i] = sm_avg;
				}
			}
			if (param->abeFnLevel > 1) {
				// 3値平均
				for (i = nSample - 1;i + 1 > 0;i--) {
					d1 = pIn[i - 0];
					d2 = pIn[i - 1];
					d3 = pIn[i - 2];
					dd1 = d1;
					dd2 = d2;
					dd3 = d3;
					dd1 >>= (56 - param->iw);
					dd2 >>= (56 - param->iw);
					dd3 >>= (56 - param->iw);
//					if (dd1 < 0) {
//						dd1 *= -1;
//					}
//					if (dd2 < 0) {
//						dd2 *= -1;
//					}
//					if (dd3 < 0) {
//						dd3 *= -1;
//					}
					dd1 /= param->abeFnLevel;
					dd2 /= param->abeFnLevel;
					dd3 /= param->abeFnLevel;
					sm_ignore = 1;
					if (dd1 + 1 == dd2 && dd2 == dd3 + 1) {
						sm_ignore = 0;
					}
					if (dd1 == dd2 + 1 && dd2 + 1 == dd3) {
						sm_ignore = 0;
					}

					if (sm_ignore == 0) {
						sm_avg = 0;
						sm_avg += pIn[i - 0];
						sm_avg += pIn[i - 1];
						sm_avg += pIn[i - 2];
						sm_avg /= 3;
						pOut[i] = sm_avg;
					}
				}
			}
			if (param->abeFnLevel > 3) {
				// 3値平均
				for (i = 0;i + 1 < nSample;i++) {
					d1 = pIn[i + 0];
					d2 = pIn[i + 1];
					d3 = pIn[i + 2];
					dd1 = d1;
					dd2 = d2;
					dd3 = d3;
					dd1 >>= (56 - param->iw);
					dd2 >>= (56 - param->iw);
					dd3 >>= (56 - param->iw);
//					if (dd1 < 0) {
//						dd1 *= -1;
//					}
//					if (dd2 < 0) {
//						dd2 *= -1;
//					}
//					if (dd3 < 0) {
//						dd3 *= -1;
//					}
					dd1 /= param->abeFnLevel;
					dd2 /= param->abeFnLevel;
					dd3 /= param->abeFnLevel;
					sm_ignore = 1;
					if (dd1 + 1 == dd2 && dd2 == dd3 + 1) {
						sm_ignore = 0;
					}
					if (dd1 == dd2 + 1 && dd2 + 1 == dd3) {
						sm_ignore = 0;
					}

					if (sm_ignore == 0) {
						sm_avg = 0;
						sm_avg += pIn[i + 0];
						sm_avg += pIn[i + 1];
						sm_avg += pIn[i + 2];
						sm_avg /= 3;
						pOut[i] = sm_avg;
					}
				}
			}
		}
#endif
#if 1
		if (param->smLowData > 0) {
			// 2値同値でその左右隣が異なる値の調整
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = 0;i + 3 < nSample;i++) {
				d1 = pIn[i + 0];
				d2 = pIn[i + 1];
				d3 = pIn[i + 2];
				d4 = pIn[i + 3];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd4 = d4;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				dd4 >>= (56 - param->iw);
				if (dd1 != dd2 && dd2 == dd3 && dd3 != dd4) {
					sm_ignore = 0;
					if ((dd1 > dd2 && (dd1 - dd2) > 2) || (dd1 < dd2 && (dd2 - dd1) > 2)) {
						sm_ignore = 1;
					}
					if ((dd3 > dd4 && (dd3 - dd4) > 2) || (dd3 < dd4 && (dd4 - dd3) > 2)) {
						sm_ignore = 1;
					}
					if (sm_ignore == 0) {
						sm_avg = (d1 + d2 + d3) / 3;
						pOut[i + 1] = sm_avg;
						sm_avg = (d2 + d3 + d4) / 3;
						pOut[i + 2] = sm_avg;
						i++;
					}
				}
			}
			// 2値同値でその左右隣が異なる値の調整(逆順)
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = nSample - 1;i > 2;i--) {
				d1 = pIn[i - 0];
				d2 = pIn[i - 1];
				d3 = pIn[i - 2];
				d4 = pIn[i - 3];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd4 = d4;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				dd4 >>= (56 - param->iw);
				if (dd1 != dd2 && dd2 == dd3 && dd3 != dd4) {
					sm_ignore = 0;
					if ((dd1 > dd2 && (dd1 - dd2) > 2) || (dd1 < dd2 && (dd2 - dd1) > 2)) {
						sm_ignore = 1;
					}
					if ((dd3 > dd4 && (dd3 - dd4) > 2) || (dd3 < dd4 && (dd4 - dd3) > 2)) {
						sm_ignore = 1;
					}
					if (sm_ignore == 0) {
						sm_avg = (d1 + d2 + d3) / 3;
						pOut[i + 1] = sm_avg;
						sm_avg = (d2 + d3 + d4) / 3;
						pOut[i + 2] = sm_avg;
						i--;
					}
				}
			}
			// 山や谷の形の波形調整
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = 0;i + 2 < nSample;i++) {
				d1 = pIn[i + 0];
				d2 = pIn[i + 1];
				d3 = pIn[i + 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				if (d1 < d2 && d2 > d3 && (dd2 - dd1) <= 2 && (dd2 - dd3) <= 2) {
					// 山
					sm_avg = (((d1 + d3) / 2) + d2) / 2;
					pOut[i + 1] = sm_avg;
				} else if (d1 > d2 && d2 < d3 && (dd1 - dd2) <= 2 && (dd3 - dd2) <= 2) {
					// 谷
					sm_avg = (((d1 + d3) / 2) + d2) / 2;
					pOut[i + 1] = sm_avg;
				}
			}
			// 山や谷の形の波形調整(逆順)
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = nSample - 1;i > 1;i--) {
				d1 = pIn[i - 0];
				d2 = pIn[i - 1];
				d3 = pIn[i - 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				if (d1 < d2 && d2 > d3 && (dd2 - dd1) <= 2 && (dd2 - dd3) <= 2) {
					// 山
					sm_avg = (((d1 + d3) / 2) + d2) / 2;
					pOut[i + 1] = sm_avg;
				} else if (d1 > d2 && d2 < d3 && (dd1 - dd2) <= 2 && (dd3 - dd2) <= 2) {
					// 谷
					sm_avg = (((d1 + d3) / 2) + d2) / 2;
					pOut[i + 1] = sm_avg;
				}
			}
			// 同値以外の移動平均
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = 0;i + 2 < nSample;i++) {
				d1 = pIn[i + 0];
				d2 = pIn[i + 1];
				d3 = pIn[i + 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				if (dd1 < 0) {
					dd1 *= -1;
				}
				if (dd2 < 0) {
					dd2 *= -1;
				}
				if (dd3 < 0) {
					dd3 *= -1;
				}
				dd1 >>= 2;
				dd2 >>= 2;
				dd3 >>= 2;
				sm_ignore = 1;
				if (dd1 == dd2 && dd2 == dd3) {
					sm_ignore = 0;
				}
				if (sm_ignore == 0) {
					sm_avg = (d1 + d2 + d3) / 3;
					pOut[i + 1] = sm_avg;
				}
			}
			// 同値以外の移動平均(逆順)
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = nSample - 1;i > 1;i--) {
				d1 = pIn[i - 0];
				d2 = pIn[i - 1];
				d3 = pIn[i - 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				if (dd1 < 0) {
					dd1 *= -1;
				}
				if (dd2 < 0) {
					dd2 *= -1;
				}
				if (dd3 < 0) {
					dd3 *= -1;
				}
				dd1 >>= 3;
				dd2 >>= 3;
				dd3 >>= 3;
				sm_ignore = 1;
				if (dd1 == dd2 && dd2 == dd3) {
					sm_ignore = 0;
				}
				if (sm_ignore == 0) {
					sm_avg = (d1 + d2 + d3) / 3;
					pOut[i + 1] = sm_avg;
				}
			}
		}
#endif
#if 1
		// 適応型フィルター処理
		// 同値が続く場合に同値の個数に応じたfftで高域カットフィルターをする。
		// フィルター後は波形がなめらかになるように制御する。
		memset(mem2,0,wkMemSize * sizeof (SSIZE));
		memset(mem3,0,wkMemSize * sizeof (SSIZE));
		memset(mem4,0,wkMemSize * sizeof (SSIZE));
		if (param->adaptiveFilder == 1) {
			pIn = mem1;							// 入力波形
			pOut = mem2;						// 同値個数
			// mem2に同値が続く場合に同値の個数を記録していく。
			for (i = 0,j = 1,cnt = 1;j < nSample;j++) {
				d1 = pIn[i] >> (56 - param->iw);
				d2 = pIn[j] >> (56 - param->iw);
				if (d1 == d2 && cnt < 255) {
					cnt++;
				} else {
					for (k = i;k < j;k++) {
						pOut[k] = cnt;
					}
					i = j;
					cnt = 1;
				}
			}
			for (i = 0;i < nSample;i++) {
				if (pOut[i] >= 0 && pOut[i] < 2) {
					pOut[i] = 0;
				}
			}
			// 同値が3つ以上続くものに対して、fft をかける
			do {
				pIn = mem1;
				pIn2 = mem2;
				memset(mem3,0,wkMemSize * sizeof (SSIZE));
				cnt = 0;
				for (i = 0;i < nSample;i++) {
					if (pIn2[i] > 0) {
						cnt = pIn2[i];
						break;
					}
				}
				if (cnt == 0) {
					break;
				}
				for (n = 0;n < 3;n++) {
					// FFT 初期設定
					copyToFFTW(fftw_in,&pIn[((fftSize / 2) * n)],fftSize);
					
					// 窓関数
					windowFFTW(fftw_in,fftSize);

					// FFT
					fftw_execute(fftw_p);

					if (param->beInfo.cutOff[cnt] > 0) {
						// 高域削除
						hfc = inSampleR / param->beInfo.cutOff[cnt];
						cutOff = ((double)fftSize / inSampleR) * hfc;
						for (i = cutOff;i < fftSize;i++) {
							fftw_out[i][0] = 0;
							fftw_out[i][1] = 0;
						}
					}
					// 半分のデータを復元
					//#pragma omp parallel for
					for (i = 1;i < fftSize / 2;i++) {
						fftw_out[fftSize - i][0] = fftw_out[i][0];
						fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
					}

					// invert FFT
					fftw_execute(fftw_ip);

					// 出力
					pOut = (SSIZE *)&mem3[((fftSize / 2) * n)];
					//#pragma omp parallel for
					for (i = 0;i < fftSize;i++) {
						pOut[i] += fftw_in[i][0] / fftSize;
					}
				}
				pIn   = (SSIZE *)mem3;
				pIn2  = (SSIZE *)mem2;
				pOut  = (SSIZE *)mem4;
				for (i = 0;i < nSample;i++) {
					if (pIn2[i] > 0 && pIn2[i] == cnt) {
						pOut[i] = pIn[i];
						pIn2[i] *= -1;
					}
				}
			} while (1);
			pIn = (SSIZE *)mem4;
			pIn2 = (SSIZE *)mem2;
			pOut = (SSIZE *)mem1;
			for (i = 0;i < nSample;i++) {
				d1 = pIn[i];
				d1 >>= (56 - param->iw);
				d2 = pOut[i];
				d2 >>= (56 - param->iw);
				if (pIn2[i] < 0) {
					if (d1 <= d2 && (d2 - d1) <= 3) {
						pOut[i] = pIn[i];
					} else if (d1 > d2 && (d1 - d2) <= 3) {
						pOut[i] = pIn[i];
					}
				}
			}
		}
#endif
		if (param->cutLowData && shresh > 0) {
			// 閾値より低いパワーの音は削除する
			pIn = mem1;
			pOut = mem3;
			memset(mem3,0,wkMemSize * sizeof (SSIZE));
			for (n = 0;n < 3;n++) {
				// FFT 初期設定
				copyToFFTW(fftw_in,&pIn[((fftSize / 2) * n)],fftSize);

				// 窓関数
				windowFFTW(fftw_in,fftSize);

				// FFT
				fftw_execute(fftw_p);

				// 削除
				for (i = 1;i < fftSize / 2;i++) {
					power = (fftw_out[i][0] * fftw_out[i][0]) + (fftw_out[i][1] * fftw_out[i][1]);
					if (power > 0) {
						power = sqrt(power);
					}
					if (power < shresh) {
						fftw_out[i][0] = 0;
						fftw_out[i][1] = 0;
					}
				}

				// 半分のデータを復元
				#pragma omp parallel for
				for (i = 1;i < fftSize / 2;i++) {
					fftw_out[fftSize - i][0] = fftw_out[i][0];
					fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
				}

				// invert FFT
				fftw_execute(fftw_ip);

				// 出力
				pOut = (SSIZE *)mem3;
				pOut = &pOut[(fftSize / 2) * n];
				#pragma omp parallel for
				for (i = 0;i < fftSize;i++) {
					pOut[i] += fftw_in[i][0] / fftSize;
				}
			}
			memcpy(mem1,mem3,wkMemSize * sizeof (SSIZE));
		}

#if 1
		// 波形の調整処理
		// 上がりや、下がりの波形の量子化誤差を少なくする
		pIn  = (SSIZE *)mem1;
		pOut = (SSIZE *)mem1;
		for (i = 4;i + 3 < nSample;) {
			next = 1;
			d1 = pIn[i - 4];
			d2 = pIn[i - 3];
			d3 = pIn[i - 2];
			d4 = pIn[i - 1];
			d5 = pIn[i];
			d6 = pIn[i + 1];
			d7 = pIn[i + 2];
			d8 = pIn[i + 3];
			if ((d2 < d3 && d3 <= d4 && d4 < d5 && d5 <= d6 && d6 < d7) ||
						(d2 <= d3 && d3 < d4 && d4 <= d5 && d5 < d6 && d6 <= d7)) {
				// 上がり波形
				samp[1] = pIn[i - 2] - pIn[i - 3];
				samp[2] = pIn[i - 1] - pIn[i - 2];
				samp[3] = pIn[i] - pIn[i - 1];
				samp[4] = pIn[i + 1] - pIn[i];
				samp[5] = pIn[i + 2] - pIn[i + 1];
				for (j = 1;j < 5;j++) {
					for (k = j + 1;k < 6;k++) {
						if (samp[j] > samp[k]) {
							samp[8] = samp[j];
							samp[j] = samp[k];
							samp[k] = samp[8];
						}
					}
				}
				samp[2] = samp[2] + samp[3] + samp[4];
				samp[2] /= 3;
				d1 = pIn[i];
				d2 = pIn[i - 1] + (SSIZE)samp[2];
				d1 >>= (56 - param->iw);
				d2 >>= (56 - param->iw);
				if (d1 == d2) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				} else if (d1 < d2 && (d2 - d1) <= 3) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				} else if (d1 > d2 && (d1 - d2) <= 3) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				}

			} else if ((d2 > d3 && d3 >= d4 && d4 > d5 && d5 >= d6 && d6 > d7) ||
						(d2 >= d3 && d3 > d4 && d4 >= d5 && d5 > d6 && d6 >= d7)) {
				// 下がり波形
				samp[1] = pIn[i - 2] - pIn[i - 3];
				samp[2] = pIn[i - 1] - pIn[i - 2];
				samp[3] = pIn[i] - pIn[i - 1];
				samp[4] = pIn[i + 1] - pIn[i];
				samp[5] = pIn[i + 2] - pIn[i + 1];
				for (j = 1;j < 5;j++) {
					for (k = j + 1;k < 6;k++) {
						if (samp[j] > samp[k]) {
							samp[8] = samp[j];
							samp[j] = samp[k];
							samp[k] = samp[8];
						}
					}
				}
				samp[2] = samp[2] + samp[3] + samp[4];
				samp[2] /= 3;
				d1 = pIn[i];
				d2 = pIn[i - 1] + (SSIZE)samp[2];
				d1 >>= (56 - param->iw);
				d2 >>= (56 - param->iw);
				if (d1 == d2) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				} else if (d1 < d2 && (d2 - d1) <= 3) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				} else if (d1 > d2 && (d1 - d2) <= 3) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				}
			}
			i += next;
		}
#endif

		if (startInSample + (fftSize / 2) >= 0) {
#if 0
			// 音のレベルに変化がないか調査
			if (param->abeNX == 1) {
				level2 = 0;
				for (i = fftSize / 2,j = 0,n = 0;n < fftSize;i++,j++,n++) {
					if (mem1[i] > 0) {
						level2 += mem1[i] >> (56 - 16);
						j++;
					}
				}
				if (j > 0) {
					level2 /= j;
				}
				if (level > 0 && level2 > 0) {
					nx = ((double)level / (double)level2);
					for (i = fftSize / 2,n = 0;n < fftSize;i++,n++) {
						mem1[i] = (SSIZE)((double)mem1[i] * nx);
					}
				}
			}
#endif
			if (outRemain >= fftSize) {
				wr = fio_write(mem1 + (fftSize / 2),sizeof (SSIZE),fftSize,fp_w);
				if (wr != fftSize) {
					param->err = STATUS_FILE_WRITE_ERR;
					return;
				}
			} else {
				wr = fio_write(mem1 + (fftSize / 2),sizeof (SSIZE),outRemain,fp_w);
				if (wr != outRemain) {
					param->err = STATUS_FILE_WRITE_ERR;
					return;
				}
			}
			if (outRemain >= fftSize) {
				outRemain -= fftSize;
			} else {
				break;
			}
		}
	}
	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	al_free(mem4);

	fftw_destroy_plan(fftw_p);
	fftw_destroy_plan(fftw_ip);
	fftw_free(fftw_in);
	fftw_free(fftw_out);
}
//---------------------------------------------------------------------------
// Function   : genNoise
// Description: 失われた高域の再現処理(正規分布のノイズ付加)
// ---
//	hfc		 	:高域のカットオフ周波数(この周波数以上の領域にデータを追加する)
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	fp_w		:出力ファイル
//	param		:変換パラメータ
//
void genNoise(long hfc,SSIZE inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3;
	long outSampleR;
	long wkMemSize;
	long fftSize,i,j,n,nn;
	
	long lowIndex,highIndex;
	long nrIndex;
	double percent,per;
	double nx;
	double p,refPw,noisePw;
	SSIZE *pIn,*pIn2,*pOut;
	SSIZE startInSample,nSample;
	
	fftw_complex *fftw_in[2],*fftw_out[2],*fftw_in2;
	fftw_plan fftw_p[2],fftw_ip[2],fftw_p2;
	double hfaNB;
	
//	// ノイズ生成用データの読み込み
//	fpNoise = fopen(param->nd_path,"rb");
//	if (fpNoise == NULL) {
//		// データがなければ作成する
//		fpNoise = fopen(param->nd_path,"wb");
//		if (fpNoise) {
//			for (i = 0;i < (long)(MAX_SAMPLE_N) * 60;i++) {
//				nrData = normalNoise() * 0x100;
//				fwrite(&nrData,1,sizeof (long),fpNoise);
//			}
//			fflush(fpNoise);
//			fclose(fpNoise);
//		}
//		fpNoise = fopen(param->nd_path,"rb");
//	}
	outSampleR = param->outSampleR;
	fio_rewind(fp_r);

	fftSize = param->outSampleR * 2;

	wkMemSize = (fftSize * 2) * sizeof (SSIZE);

	mem1 = (SSIZE *)al_malloc(wkMemSize);
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem2 = (SSIZE *)al_malloc(wkMemSize);
	if (mem2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem3 = (SSIZE *)al_malloc(wkMemSize);
	if (mem3 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_in[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize / 100);
	if (fftw_in2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_p[0] = fftw_plan_dft_1d(fftSize,fftw_in[0],fftw_out[0],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p[1] = fftw_plan_dft_1d(fftSize,fftw_in[1],fftw_out[1],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[0] = fftw_plan_dft_1d(fftSize,fftw_out[0],fftw_in[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[1] = fftw_plan_dft_1d(fftSize,fftw_out[1],fftw_in[1],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p2 = fftw_plan_dft_1d(fftSize / 100,fftw_in2,fftw_in2,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	per = -1;
	nrIndex = 0;
	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + (fftSize + (fftSize / 2));startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
			percent = ((double)startInSample / inSample);
			percent *= 100;
			percent = (double)((int)percent);
			if (percent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)percent);
				fflush(stdout);
				if (chkAbort(param,percent,10) == 1) exit(0);
			}
			per = percent;
		}

		memset(mem1,0,wkMemSize);
		memset(mem3,0,wkMemSize);

		// mem2 には正規分布のノイズを格納する
		pIn2 = (SSIZE *)mem2;
		#pragma omp parallel for
		for (i = 0;i < fftSize * 2;i++) {
			pIn2[i] = normalNoise() * 0x100;
		}

		if (startInSample >= 0 && startInSample + (fftSize * 2) < inSample) {
			nSample = fftSize * 2;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 2;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}

		pIn = (SSIZE *)mem1;
		for (n = 0;n < 3;n++) {
			lowIndex = ((double)fftSize / outSampleR) * (hfc - 2000);
			highIndex = ((double)fftSize / outSampleR) * hfc;
			// FFT 初期設定
			copyToFFTW(fftw_in[0],&pIn[((fftSize / 2) * n)],fftSize);
			windowFFTW(fftw_in[0],fftSize);

			// FFT
			fftw_execute(fftw_p[0]);

			// 元信号の高域のパワーを調べる
			refPw = 0;
			for (i = lowIndex,nn = 0;i < highIndex;i++,nn++) {
				p = fftw_out[0][i][0] * fftw_out[0][i][0] + fftw_out[0][i][1] * fftw_out[0][i][1];
				if (p != 0) {
					p = sqrt(p);
				}
				refPw += p;
			}
			if (nn > 0) {
				refPw /= nn;
			}

			// 付加する信号
			copyToFFTW(fftw_in[1],&pIn2[((fftSize / 2) * n)],fftSize);
			windowFFTW(fftw_in[1],fftSize);

			// FFT
			fftw_execute(fftw_p[1]);
			// 1/f 信号にする
			adjPinkFilter(1,fftSize,fftw_out[1],param);

			// 付加する信号のパワーを調べる
			noisePw = 0;
			lowIndex = ((double)fftSize / outSampleR) * hfc;
			highIndex = ((double)fftSize / outSampleR) * (hfc + 2000);
			for (i = lowIndex,nn = 0;i < highIndex;i++,nn++) {
				p = fftw_out[1][i][0] * fftw_out[1][i][0] + fftw_out[1][i][1] * fftw_out[1][i][1];
				if (p != 0) {
					p = sqrt(p);
				}
				noisePw += p;
			}
			if (nn > 0) {
				noisePw /= nn;
			}
			if (noisePw == 0) {
				noisePw = 1;
			}
			nx = refPw / noisePw;
			highIndex = ((double)fftSize / outSampleR) * hfc;

			#pragma omp parallel for
			for (i = highIndex;i < fftSize / 2;i++) {
				fftw_out[1][i][0] *= nx;
				fftw_out[1][i][1] *= nx;
			}

			if (param->hfa != 1 && param->hfaNB > 0) {
				hfaNB = param->hfaNB / 100.0;
				// hfa1 信号
				for (i = highIndex;i < fftSize / 2;i++) {
					fftw_out[1][i][0] *= hfaNB;
					fftw_out[1][i][1] *= hfaNB;
				}
				hfaNB = 1.0 - hfaNB;
				// hfa2 信号
				for (i = highIndex;i < fftSize / 2;i++) {
					fftw_out[0][i][0] *= hfaNB;
					fftw_out[0][i][1] *= hfaNB;
				}
				// 合成
				for (i = highIndex;i < fftSize / 2;i++) {
					fftw_out[1][i][0] += fftw_out[0][i][0];
					fftw_out[1][i][1] += fftw_out[0][i][1];
				}
			}

			// 低域カット
			highIndex = ((double)fftSize / outSampleR) * hfc;

			#pragma omp parallel for
			for (i = 1;i < highIndex;i++) {
				fftw_out[1][i][0] = 0;
				fftw_out[1][i][1] = 0;
			}

			// 半分のデータを復元
			#pragma omp parallel for
			for (i = 1;i < fftSize / 2;i++) {
				fftw_out[1][fftSize - i][0] = fftw_out[1][i][0];
				fftw_out[1][fftSize - i][1] = fftw_out[1][i][1] * -1;
			}
			// 直流部分除去
			fftw_out[1][0][0] = 0;
			fftw_out[1][0][1] = 0;

			// invert FFT
			fftw_execute(fftw_ip[1]);

			pOut = (SSIZE *)mem3;
			pOut = &pOut[(fftSize / 2) * n];
			#pragma omp parallel for
			for (i = 0;i < fftSize;i++) {
				pOut[i] += fftw_in[1][i][0] / fftSize;
			}
		}
		if (startInSample + fftSize / 2 >= 0) {
			if (param->hfa == 1) {
				pIn = (SSIZE *)mem1;
			} else {
				pIn = (SSIZE *)mem3;
			}
			pIn = pIn + (fftSize / 2);
			pOut = (SSIZE *)mem3;
			pOut = pOut + (fftSize / 2);
			for (i = 0;i < 100;i++) {
				// 元信号に音があるかを調べる
				for (j = 0;j < fftSize / 100;j++) {
					fftw_in2[j][0] = pIn[j];
					fftw_in2[j][1] = 0;
				}
				// 窓関数
				for (j = 0;j < ((fftSize / 100) - 1) / 2;j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 * j / ((double)fftSize / 100));
				}
				for (j = ((fftSize / 100) - 1) / 2;j < (fftSize / 100);j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 - 2.0 * j / ((double)fftSize / 100));
				}

				// FFT
				fftw_execute(fftw_p2);

				// 元信号の高域のパワーを調べる
				refPw = 0;
				lowIndex = (((double)fftSize / 100) / outSampleR) * (hfc - 2000);
				highIndex = (((double)fftSize / 100) / outSampleR) * hfc;
				for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
					p = fftw_in2[j][0] * fftw_in2[j][0] + fftw_in2[j][1] * fftw_in2[j][1];
					if (p != 0) {
						p = sqrt(p);
					}
					refPw += p;
				}
				if (nn > 0) {
					refPw /= nn;
				}
				// 付加信号の高域のパワーを調べる
				for (j = 0;j < fftSize / 100;j++) {
					fftw_in2[j][0] = pOut[j];
					fftw_in2[j][1] = 0;
				}
				// 窓関数
				for (j = 0;j < ((fftSize / 100) - 1) / 2;j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 * j / ((double)fftSize / 100));
				}
				for (j = ((fftSize / 100) - 1) / 2;j < (fftSize / 100);j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 - 2.0 * j / ((double)fftSize / 100));
				}

				// FFT
				fftw_execute(fftw_p2);

				noisePw = 0;
				lowIndex = (((double)fftSize / 100) / outSampleR) * hfc;
				highIndex = (((double)fftSize / 100) / outSampleR) * (hfc + 2000);
				for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
					p = fftw_in2[j][0] * fftw_in2[j][0] + fftw_in2[j][1] * fftw_in2[j][1];
					if (p != 0) {
						p = sqrt(p);
					}
					noisePw += p;
				}
				if (nn > 0) {
					noisePw /= nn;
				}
				if (refPw > 0) {
					if (param->hfa == 1) {
						for (j = 0;j < fftSize / 100;j++) {
							pOut[j] *= ((refPw / noisePw) * 0.70);
						}
					} else {
						for (j = 0;j < fftSize / 100;j++) {
							pOut[j] *= (refPw / noisePw);
						}
					}
				} else {
					for (j = 0;j < fftSize / 100;j++) {
						pOut[j] = 0;
					}
				}
				pIn  += fftSize / 100;
				pOut += fftSize / 100;
			}
			// 再び低域カット処理をする
			pIn = (SSIZE *)mem1;
			pOut = (SSIZE *)mem3;
			for (n = 0;n < 3;n++) {
				// FFT 初期設定
				for (i = 0;i < fftSize;i++) {
					fftw_in[1][i][0] = pOut[((fftSize / 2) * n) + i];
				   	fftw_in[1][i][1] = 0;
				}
				// 窓関数
				for (i = 0;i < (fftSize - 1) / 2;i++) {
					fftw_in[1][i][0] = fftw_in[1][i][0] * (2.0 * i / (double)fftSize);
				}
				for (i = (fftSize - 1) / 2;i < fftSize;i++) {
					fftw_in[1][i][0] = fftw_in[1][i][0] * (2.0 - 2.0 * i / (double)fftSize);
				}

				// FFT
				fftw_execute(fftw_p[1]);

				highIndex = ((double)fftSize / outSampleR) * (hfc - 1000);

				// 低域カット
				for (i = 1;i < highIndex;i++) {
					fftw_out[1][i][0] = 0;
					fftw_out[1][i][1] = 0;
				}

				// 半分のデータを復元
				for (i = 1;i < fftSize / 2;i++) {
					fftw_out[1][fftSize - i][0] = fftw_out[1][i][0];
					fftw_out[1][fftSize - i][1] = fftw_out[1][i][1] * -1;
				}

				// 直流部分除去
				fftw_out[1][0][0] = 0;
				fftw_out[1][0][1] = 0;

				// invert FFT
				fftw_execute(fftw_ip[1]);

				// 出力
				for (i = 0;i < fftSize;i++) {
					pOut[((fftSize / 2) * n) + i] += fftw_in[1][i][0] / fftSize;
				}
			}
			for (i = 0;i < fftSize;i++) {
				pOut[(fftSize / 2) + i] += pIn[(fftSize / 2) + i];
			}
			fio_seek(fp_w,(startInSample + (fftSize / 2)) * sizeof (SSIZE),SEEK_SET);
			outTempFile(fp_w,mem3 + fftSize / 2,fftSize,param);
			if (param->err) {
				break;
			}
		}
	}

	fio_flush(fp_w);

	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	fftw_destroy_plan(fftw_p[0]);
	fftw_destroy_plan(fftw_p[1]);
	fftw_destroy_plan(fftw_ip[0]);
	fftw_destroy_plan(fftw_ip[1]);
	fftw_destroy_plan(fftw_p2);
	fftw_free(fftw_in[0]);
	fftw_free(fftw_in[1]);
	fftw_free(fftw_out[0]);
	fftw_free(fftw_out[1]);
	fftw_free(fftw_in2);

}
//---------------------------------------------------------------------------
// Function   : genOverTone
// Description: 失われた高域の再現処理(倍音解析)
// ---
//	hfc		 	:高域のカットオフ周波数(この周波数以上の領域にデータを追加する)
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	fp_w		:出力ファイル
//	param		:変換パラメータ
//
void genOverTone(long hfc,SSIZE inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param,FFT_PARAM *fft_param)
{
	long outSampleR;
	SSIZE *mem0,*mem1,*mem2,*mem3,*mem4;
	long wkMemSize;
	long fftSize,i,j,k,nn;
	long lowIndex,highIndex;
	double percent,per;
	double nx;
	double p,refPw,ovTonePw;
	double avg;
	SSIZE *pIn[4],*pOut[4];
	SSIZE startInSample,nSample;
	SSIZE s[2];
	SSIZE p1,p2,p3;
	fftw_complex *fftw_in[3],*fftw_out[3],*fftw_in2;
	fftw_plan fftw_p[3],fftw_ip[3],fftw_p2;
	OVERTONE_INFO *ovInfo[3];
	double *phaseX,*phaseY;
	outSampleR = param->outSampleR;
	fio_rewind(fp_r);
	fio_rewind(fp_w);
param->hfa3_max = 0;

	fftSize = 4096;
	fftSize = outSampleR / 14;
	if (outSampleR == 44100 * 2 || outSampleR == 48000 * 2) {
		fftSize = outSampleR / 10;
		fftSize = 4096 * 2;
		fftSize = outSampleR / 14;
	}
	if (outSampleR == 32000 * 6 || outSampleR == 44100 * 4 || outSampleR == 48000 * 4) {
		fftSize = outSampleR / 10;
		fftSize = 8192 * 2;
		fftSize = outSampleR / 14;
	}
	if (outSampleR == 32000 * 12 || outSampleR == 44100 * 8 || outSampleR == 48000 * 8) {
		fftSize = outSampleR / 10;
		fftSize = 16384 * 2;
		fftSize = outSampleR / 16;
	}
	if (outSampleR == 32000 * 24 || outSampleR == 44100 * 16 || outSampleR == 48000 * 16) {
		fftSize = 32768 * 2;
		fftSize = outSampleR / 16;
	}
	if (outSampleR == 32000 * 48 || outSampleR == 44100 * 32 || outSampleR == 48000 * 32) {
		fftSize = 65536 * 2;
		fftSize = outSampleR / 18;
	}
	if (outSampleR == 32000 * 64 || outSampleR == 44100 * 64 || outSampleR == 48000 * 64) {
		fftSize = 65536 * 2;
		fftSize = outSampleR / 20;
	}

//fftSize = 4096;
	wkMemSize = (fftSize * 2) * sizeof (SSIZE);

	mem1 = (SSIZE *)al_malloc(wkMemSize);
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem2 = (SSIZE *)al_malloc(wkMemSize);
	if (mem2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem3 = (SSIZE *)al_malloc(wkMemSize);
	if (mem3 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem4 = (SSIZE *)al_malloc(wkMemSize);
	if (mem4 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	phaseX = (double *)al_malloc(65536 * 2 * sizeof (double));
	if (phaseX == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	phaseY = (double *)al_malloc(65536 * 2 * sizeof (double));
	if (phaseY == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	ovInfo[0] = (OVERTONE_INFO *)malloc(sizeof (OVERTONE_INFO));
	if (ovInfo[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	ovInfo[1] = (OVERTONE_INFO *)malloc(sizeof (OVERTONE_INFO));
	if (ovInfo[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	ovInfo[2] = (OVERTONE_INFO *)malloc(sizeof (OVERTONE_INFO));
	if (ovInfo[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 1
	fftw_in[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_in[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_out[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 2
	fftw_in[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_in[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_out[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 3
	fftw_in[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_in[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_out[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 1
	fftw_p[0] = fftw_plan_dft_1d(fftSize,fftw_in[0],fftw_out[0],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[0] = fftw_plan_dft_1d(fftSize,fftw_out[0],fftw_in[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 2
	fftw_p[1] = fftw_plan_dft_1d(fftSize,fftw_in[1],fftw_out[1],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[1] = fftw_plan_dft_1d(fftSize,fftw_out[1],fftw_in[1],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 3
	fftw_p[2] = fftw_plan_dft_1d(fftSize,fftw_in[2],fftw_out[2],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[2] = fftw_plan_dft_1d(fftSize,fftw_out[2],fftw_in[2],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_in2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p2 = fftw_plan_dft_1d(fftSize/12,fftw_in2,fftw_in2,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	avg = 0;
	for (i = 1;i < fftSize / 2;i++) {
		double t = ((2 * M_PI) / ((fftSize / 2) - 1)) * i;
		phaseX[i] = sin(t) * 1000000;
		phaseY[i] = cos(t) * 1000000;
		p = (phaseX[i] * phaseX[i]) + (phaseY[i] * phaseY[i]);
		if (p != 0) {
			p = sqrt(p);
		}
		avg += p;
	}
	if (avg != 0) {
		avg /= (fftSize / 2);
	}
	for (i = 1;i < fftSize / 2;i++) {
		p = (phaseX[i] * phaseX[i]) + (phaseY[i] * phaseY[i]);
		if (p != 0) {
			p = sqrt(p);
			nx = avg / p;
			phaseX[i] *= nx;
			phaseY[i] *= nx;
		}
	}

	memset(ovInfo[0],0,sizeof (OVERTONE_INFO));
	memset(ovInfo[1],0,sizeof (OVERTONE_INFO));
	memset(ovInfo[2],0,sizeof (OVERTONE_INFO));

	for (i = 0;i < 360;i++) {
		ovInfo[0]->do2Idx[i] = -1;
		ovInfo[1]->do2Idx[i] = -1;
		ovInfo[2]->do2Idx[i] = -1;
	}
	ovInfo[0]->nSample = fftSize / 2;
	ovInfo[0]->validSamplingRate = hfc;
	ovInfo[0]->samplingRate = outSampleR;
	ovInfo[0]->phaseX = phaseX;
	ovInfo[0]->phaseY = phaseY;

	ovInfo[1]->nSample = fftSize / 2;
	ovInfo[1]->validSamplingRate = hfc;
	ovInfo[1]->samplingRate = outSampleR;
	ovInfo[1]->phaseX = phaseX;
	ovInfo[1]->phaseY = phaseY;
	ovInfo[1]->log = 1;

	ovInfo[2]->nSample = fftSize / 2;
	ovInfo[2]->validSamplingRate = hfc;
	ovInfo[2]->samplingRate = outSampleR;
	ovInfo[2]->phaseX = phaseX;
	ovInfo[2]->phaseY = phaseY;

	pIn[0]	= &mem1[((fftSize / 2) * 0)];
	pOut[0] = &mem2[((fftSize / 2) * 0)];
	pIn[1]	= &mem1[((fftSize / 2) * 1)];
	pOut[1] = &mem3[((fftSize / 2) * 1)];
	pIn[2]	= &mem1[((fftSize / 2) * 2)];
	pOut[2] = &mem4[((fftSize / 2) * 2)];

	per = -1;
	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + ((fftSize / 2));startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
			percent = ((double)startInSample / inSample);
			percent *= 100;
			percent = (double)((int)percent);
			if (percent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)percent);
				fflush(stdout);
				if (chkAbort(param,percent,3) == 1) exit(0);
			}
			per = percent;
		}
		memset(mem1,0,wkMemSize);

		if (startInSample >= 0 && startInSample + (fftSize * 2) < inSample) {
			nSample = fftSize * 2;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 2;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}

		memset(mem2,0,wkMemSize);
		memset(mem3,0,wkMemSize);
		memset(mem4,0,wkMemSize);
		
		#pragma omp parallel
		{
			#pragma omp sections
			{
				#pragma omp section
				{
					// 1
					genOverToneSub(hfc,pIn[0],pOut[0],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],ovInfo[0],param,fft_param);
				}
				#pragma omp section
				{
					// 2
					genOverToneSub(hfc,pIn[1],pOut[1],fftw_in[1],fftw_out[1],fftw_p[1],fftw_ip[1],ovInfo[1],param,fft_param);
				}
				#pragma omp section
				{
					// 3
					genOverToneSub(hfc,pIn[2],pOut[2],fftw_in[2],fftw_out[2],fftw_p[2],fftw_ip[2],ovInfo[2],param,fft_param);
				}
			}
		}

		for (i = 0;i < fftSize * 2;i++) {
			mem2[i] += mem3[i] + mem4[i];
		}

		if (startInSample + fftSize / 2 >= 0) {
			SSIZE *pIn2,*pOut2;
			// レベル調整
			pIn2  = (SSIZE *)mem1;
			pOut2 = (SSIZE *)mem2;
			pIn2  += (fftSize / 2);
			pOut2 += (fftSize / 2);
			for (i = 0;i < 12;i++) {
				for (j = 0;j < fftSize / 12;j++) {
					fftw_in2[j][0] = pIn2[j];
					fftw_in2[j][1] = 0;
				}
				pIn2 += fftSize / 12;
				// 窓関数
				for (j = 0;j < ((fftSize / 12) - 1) / 2;j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 * j / ((double)fftSize / 12));
				}
				for (j = ((fftSize / 12) - 1) / 2;j < (fftSize / 12);j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 - 2.0 * j / ((double)fftSize / 12));
				}

				// FFT
				fftw_execute(fftw_p2);

				// 元信号の高域のパワーを調べる
				refPw = 0;
				lowIndex = (((double)fftSize / 12) / outSampleR) * (hfc - 2000);
				highIndex = (((double)fftSize / 12) / outSampleR) * hfc;
				memset(ovInfo[0]->avg,0,sizeof (double));
				for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
					p = fftw_in2[j][0] * fftw_in2[j][0] + fftw_in2[j][1] * fftw_in2[j][1];
					if (p != 0) {
						p = sqrt(p);
					}
					ovInfo[0]->avg[nn] = p;
				}
				for (j = 0;j + 1 < nn;j++) {
					for (k = j;k < nn;k++) {
						if (ovInfo[0]->avg[j] > ovInfo[0]->avg[k]) {
							p = ovInfo[0]->avg[j];
							ovInfo[0]->avg[j] = ovInfo[0]->avg[k];
							ovInfo[0]->avg[k] = p;
						}
					}
				}
				for (j = 0;j + 2 < nn;j++) {
					refPw += ovInfo[0]->avg[j];
				}
				if (j > 0) {
					refPw /= j;
				}

				// 付加信号の高域のパワーを調べる
				for (j = 0;j < fftSize / 12;j++) {
					fftw_in2[j][0] = pOut2[j];
					fftw_in2[j][1] = 0;
				}
				// 窓関数
				for (j = 0;j < ((fftSize / 12) - 1) / 2;j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 * j / ((double)fftSize / 12));
				}
				for (j = ((fftSize / 12) - 1) / 2;j < (fftSize / 12);j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 - 2.0 * j / ((double)fftSize / 12));
				}
				// FFT
				fftw_execute(fftw_p2);

				ovTonePw = 0;
				lowIndex = (((double)fftSize / 12) / outSampleR) * hfc;
				highIndex = (((double)fftSize / 12) / outSampleR) * (hfc + 2000);
				memset(ovInfo[0]->avg,0,sizeof (double));
				for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
					p = fftw_in2[j][0] * fftw_in2[j][0] + fftw_in2[j][1] * fftw_in2[j][1];
					if (p != 0) {
						p = sqrt(p);
					}
					ovInfo[0]->avg[nn] = p;
				}
				for (j = 0;j + 1 < nn;j++) {
					for (k = j;k < nn;k++) {
						if (ovInfo[0]->avg[j] > ovInfo[0]->avg[k]) {
							p = ovInfo[0]->avg[j];
							ovInfo[0]->avg[j] = ovInfo[0]->avg[k];
							ovInfo[0]->avg[k] = p;
						}
					}
				}
				for (j = 0;j + 3 < nn;j++) {
					ovTonePw += ovInfo[0]->avg[j];
				}
				if (j > 0) {
					ovTonePw /= j;
				}
				if (ovTonePw > 0) {
					for (j = 0;j < fftSize / 12;j++) {
						pOut2[j] *= ((refPw / ovTonePw));
						pOut2[j] *= 0.95;
					}
				} else {
					for (j = 0;j < fftSize / 12;j++) {
						pOut2[j] = 0;
					}
				}
				pOut2 += fftSize / 12;
			}

#if 1
			// パルス除去
			pIn2  = (SSIZE *)mem2;
			for (i = 1;i + 1 < fftSize;i++) {
				p1 = pIn2[i - 1];
				p2 = pIn2[i];
				p3 = pIn2[i + 1];
				p1 >>= (56 - param->w);
				p2 >>= (56 - param->w);
				p3 >>= (56 - param->w);
				if (p1 < p2 && p2 > p3) {
					if (p2 - ((p1 + p3) / 2) >= 4) {
						pIn2[i] = ((pIn2[i - 1] + pIn2[i + 1]) / 2);
					}
				} else if (p1 > p2 && p2 < p3) {
					if (((p1 + p3) / 2) - p2 >= 4) {
						pIn2[i] = ((pIn2[i - 1] + pIn2[i + 1]) / 2);
					}
				}
			}
#endif
			outTempFile(fp_w,mem2 + fftSize / 2,fftSize,param);
			if (param->err) {
				break;
			}
		}
	}

	fio_flush(fp_w);
	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	al_free(mem4);

	al_free(phaseX);
	al_free(phaseY);

	free(ovInfo[0]);
	free(ovInfo[1]);
	free(ovInfo[2]);

	fftw_destroy_plan(fftw_p[0]);
	fftw_destroy_plan(fftw_ip[0]);

	fftw_destroy_plan(fftw_p[1]);
	fftw_destroy_plan(fftw_ip[1]);

	fftw_destroy_plan(fftw_p[2]);
	fftw_destroy_plan(fftw_ip[2]);

	fftw_destroy_plan(fftw_p2);

	fftw_free(fftw_in[0]);
	fftw_free(fftw_out[0]);

	fftw_free(fftw_in[1]);
	fftw_free(fftw_out[1]);

	fftw_free(fftw_in[2]);
	fftw_free(fftw_out[2]);

	fftw_free(fftw_in2);
}

//---------------------------------------------------------------------------
// Function   : genOverToneSub
// Description: 失われた高域の再現処理のサブ関数(倍音解析)
// ---
//	hfc		 	:高域のカットオフ周波数(この周波数以上の領域にデータを追加する)
//	pIn			:入力バッファ
//	pOut		:出力バッファ
//	fftw_in		:FFTW 入力
//	fftw_out	:FFTW 出力
//	fftw_p		:FFTW プラン
//	fftw_ip		:FFTW プラン
//	ovInfo		:高域生成用構造体
//	param		:変換パラメータ
//
void genOverToneSub(long hfc,SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,OVERTONE_INFO *ovInfo,PARAM_INFO *param,FFT_PARAM *fft_param)
{
	long outSampleR;
	long fftSize,i,j,n,nn;
	long lowIndex,highIndex;
	double nx;
	double p,p2,refPw,ovTonePw;
	double phaseTemp;
	int d,d2;
	int overToneNotFound;

	outSampleR = param->outSampleR;
	fftSize = outSampleR / 14;
	if (outSampleR == 44100 * 2 || outSampleR == 48000 * 2) {
		// 88200,96000
		fftSize = outSampleR / 10;
		fftSize = 4096 * 2;
		fftSize = outSampleR / 14;
	}
	if (outSampleR == 32000 * 6 || outSampleR == 44100 * 4 || outSampleR == 48000 * 4) {
		// 128000,176400,192000
		fftSize = outSampleR / 10;
		fftSize = 8192 * 2;
		fftSize = outSampleR / 14;
	}
	if (outSampleR == 32000 * 12 || outSampleR == 44100 * 8 || outSampleR == 48000 * 8) {
		// 352800,384000
		fftSize = outSampleR / 10;
		fftSize = 16384 * 2;
		fftSize = outSampleR / 16;
	}
	if (outSampleR == 32000 * 24 || outSampleR == 44100 * 16 || outSampleR == 48000 * 16) {
		// 705600,768000
		fftSize = 32768 * 2;
		fftSize = outSampleR / 16;
	}
	if (outSampleR == 32000 * 48 || outSampleR == 44100 * 32 || outSampleR == 48000 * 32 || outSampleR == 44100 * 64 || outSampleR == 48000 * 64) {
		// 1411200,1536000
		fftSize = 65536 * 2;
		fftSize = outSampleR / 18;
	}
	if (outSampleR == 32000 * 64 || outSampleR == 44100 * 64 || outSampleR == 48000 * 64) {
		// 
		fftSize = 65536 * 2;
		fftSize = outSampleR / 20;
	}


//fftSize = 4096;
	// FFT 初期設定
	copyToFFTW(fftw_in,pIn,fftSize);

	// 窓関数
	windowFFTW(fftw_in,fftSize);

	// FFT
	fftw_execute(fftw_p);

	// 元信号の高域のパワーを調べる
	memset(ovInfo->pw_cnt,0,65536*2 * sizeof (int));
	refPw = 0;
	lowIndex = ((double)fftSize / outSampleR) * (hfc - 2000);
	highIndex = ((double)fftSize / outSampleR) * hfc;
#if 0
	for (i = lowIndex,nn = 0;i < highIndex;i++,nn++) {
		p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
		if (p != 0) {
			p = sqrt(p);
		}
		refPw += p;
	}
#endif
	for (i = lowIndex;i < highIndex;i++) {
		p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
		if (p != 0) {
			p = sqrt(p);
			p = log10(p) * 20;
			p /= 3;
		}
		for (j = lowIndex;j < highIndex;j++) {
			p2 = fftw_out[j][0] * fftw_out[j][0] + fftw_out[j][1] * fftw_out[j][1];
			if (p2 != 0) {
				p2 = sqrt(p2);
				p2 = log10(p2) * 20;
				p2 /= 3;
			}
			if (p > 0 || p2 > 0) {
				if ((long)p == (long)p2) {
					ovInfo->pw_cnt[i]++;
				}
			}
		}
	}
	for (i = lowIndex + 1,nn = lowIndex;i < highIndex;i++) {
		if (ovInfo->pw_cnt[i] > ovInfo->pw_cnt[nn]) {
			nn = i;
		}
	}
	refPw = fftw_out[nn][0] * fftw_out[nn][0] + fftw_out[nn][1] * fftw_out[nn][1];
	if (refPw > 0) {
		refPw = sqrt(refPw);
	}

	//
	// 倍音解析
	memset(ovInfo->power,0,65536*2 * sizeof (double));
	memset(ovInfo->phase,0,65536*2 * sizeof (double));
	memset(ovInfo->pw,0,65536*2 * sizeof (double));
	memset(ovInfo->avg,0,65536*2 * sizeof (double));
	memset(ovInfo->diff,0,65536*2 * sizeof (double));
	memset(ovInfo->base,0,65536*2 * sizeof (double));
	memset(ovInfo->baseToAdj,0,65536*2 * sizeof (double));

	for (i = 1;i < fftSize / 2;i++) {
		p = (fftw_out[i][0] * fftw_out[i][0]) + (fftw_out[i][1] * fftw_out[i][1]);
		if (p != 0) {
			p = sqrt(p);
		}
		if (fftw_out[i][0] != 0 && fftw_out[i][1] != 0) {
			phaseTemp = atan2(fftw_out[i][1],fftw_out[i][0]) * 180 / M_PI;
			phaseTemp += 180;
		} else {
			phaseTemp = 0;
		}
		ovInfo->phase[i] = phaseTemp;
		ovInfo->power[i] = p;
	}
#if 0
if (param->r1_flag && ovInfo->log) {
	FILE *fp;
	fp = fopen("d:\\fft_p.csv","a");
	if (fp) {
		fprintf(fp,"\n");
		for (i = 1;i < fftSize / 2;i++) {
			fprintf(fp,",%lf",ovInfo->power[i]);
		}
		fclose(fp);
	}
	fp = fopen("d:\\fft_h.csv","a");
	if (fp) {
		fprintf(fp,"\n");
		for (i = 1;i < fftSize / 2;i++) {
			fprintf(fp,",%lf",ovInfo->phase[i]);
		}
		fclose(fp);
	}
}
#endif

#if 0
if (ovInfo->log) {
	FILE *ofp;
	ofp = fopen("d:\\fft.csv","a");
	if (ofp) {
		for (i = 300;i < 750;i++) {
			fprintf(ofp,"%lf,",ovInfo->power[i]);
		}
		fprintf(ofp,"\n");
		fclose(ofp);
	}
}
#endif
	if (param->hfa == 2) {
		anaOverToneHFA2(ovInfo,param);
	} else {
		anaOverToneHFA3(ovInfo,param);
	}

	//
	// 信号生成
	for (i = 1;i < fftSize / 2;i++) {
		if (ovInfo->pw[i] > 0) {
			d = ovInfo->phase[i];
			if (d < 0 || d >= 360) {
				d = 0;
			}
			if (ovInfo->do2Idx[d] == -1) {
				for (j = 1;j < fftSize / 2;j++) {
					if (ovInfo->phaseX[j] != 0 && ovInfo->phaseY[j] != 0) {
						phaseTemp = atan2(ovInfo->phaseY[j],ovInfo->phaseX[j]) * 180 / M_PI;
						phaseTemp += 180;
					} else {
						phaseTemp = 0;
					}
					d2 = phaseTemp;
					if (d == d2) {
						break;
					}
				}
			} else {
				j = ovInfo->do2Idx[d];
			}
			if (j < fftSize / 2) {
				// 位相が一致
				fftw_out[i][0] = ovInfo->phaseX[j];
				fftw_out[i][1] = ovInfo->phaseY[j];
				ovInfo->do2Idx[d] = j;
			} else {
				ovInfo->pw[i] = 0;
			}
		}
	}
	overToneNotFound = 1;
	for (i = 1;i < fftSize / 2;i++) {
		if (ovInfo->pw[i] > 0) {
			fftw_out[i][0] *= ovInfo->pw[i];
			fftw_out[i][1] *= ovInfo->pw[i];
			overToneNotFound = 0;
		} else {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
	}

	memset(ovInfo->pw_cnt,0,65536*2 * sizeof (int));
	lowIndex = ((double)fftSize / outSampleR) * (hfc);
	highIndex = ((double)fftSize / outSampleR) * (hfc + 2000);
	for (i = lowIndex;i < highIndex;i++) {
		p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
		if (p != 0) {
			p = sqrt(p);
			p = log10(p) * 20;
			p /= 3;
		}

		for (j = lowIndex;j < highIndex;j++) {
			p2 = fftw_out[j][0] * fftw_out[j][0] + fftw_out[j][1] * fftw_out[j][1];
			if (p2 != 0) {
				p2 = sqrt(p2);
				p2 = log10(p2) * 20;
				p2 /= 3;
			}
			if (p > 0 || p2 > 0) {
				if ((long)p == (long)p2) {
					ovInfo->pw_cnt[i]++;
				}
			}
		}
	}
	for (i = lowIndex + 1,nn = lowIndex;i < highIndex;i++) {
		if (ovInfo->pw_cnt[i] > ovInfo->pw_cnt[nn]) {
			nn = i;
		}
	}
	ovTonePw = fftw_out[nn][0] * fftw_out[nn][0] + fftw_out[nn][1] * fftw_out[nn][1];
	if (ovTonePw > 0) {
		ovTonePw = sqrt(ovTonePw);
	} else {
		ovTonePw = 1;
	}
	if (overToneNotFound == 0) {
		// 付加する信号のパワーを調べる
#if 0
		ovTonePw = 0;
		for (i = lowIndex,nn = 0;i < highIndex;i++,nn++) {
			p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
			if (p != 0) {
				p = sqrt(p);
			}
			ovTonePw += p;
		}
		if (nn > 0) {
			ovTonePw /= nn;
		}
		if (ovTonePw == 0) {
			ovTonePw = 1;
		}
#endif
		nx = refPw / ovTonePw;
		nx *= 0.9;

		for (i = 1;i < fftSize / 2;i++) {
			fftw_out[i][0] *= nx;
			fftw_out[i][1] *= nx;
		}
		adjPinkFilter(2,fftSize,fftw_out,param);

		lowIndex = ((double)fftSize / outSampleR) * (hfc - 2000);
		highIndex = ((double)fftSize / outSampleR) * hfc;

		// 低域カット
		for (i = 1;i <= highIndex;i++) {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
		// 半分のデータを復元
		for (i = 1;i < fftSize / 2;i++) {
			fftw_out[fftSize - i][0] = fftw_out[i][0];
			fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
		}

		// 直流成分を削除
		fftw_out[0][0] = 0;
		fftw_out[0][1] = 0;

		// invert FFT
		fftw_execute(fftw_ip);

		// 窓関数
		for (i = 0;i < (fftSize - 1) / 2;i++) {
			fftw_in[i][0] = fftw_in[i][0] * (2.0 * i / (double)fftSize);
		}
		for (i = (fftSize - 1) / 2;i < fftSize;i++) {
			fftw_in[i][0] = fftw_in[i][0] * (2.0 - 2.0 * i / (double)fftSize);
		}
	} else {
		// 無音生成
		for (i = 0;i < fftSize;i++) {
			fftw_in[i][0] = 0;
		}
	}
	// 出力
	for (i = 0;i < fftSize;i++) {
		pOut[i] += (fftw_in[i][0] / fftSize);
	}
}

//---------------------------------------------------------------------------
// Function   : anaOverToneHFA2
// Description: 倍音解析
//				window幅ごとにpowerの平均値をとりpowerが大きいものを抽出する
// ---
//	ovInfo :倍音生成用構造体
//	param  :パラメータ
//
void anaOverToneHFA2(OVERTONE_INFO *ovInfo,PARAM_INFO *param)
{
	DWORD ofs,window,width,swidth;
	DWORD n,i;
	double avg,maxAvg;
	double avgLine;
	int maxOfs,maxWin;
	
	
	int lowIndex,highIndex;
	int pha;
	
	
	long nSample;
	long lowHz,wid;
	long skipCnt;
	nSample = ovInfo->nSample;
	//
	// 計算対象のインデックスを求める
	for (i = 1;i < ovInfo->nSample;i++) {
		ovInfo->pw[i] = 0;
	}

	lowHz = 9000;
	wid   = 2500;
	if (lowHz + wid > ovInfo->validSamplingRate) {
		lowHz = ovInfo->validSamplingRate - wid;
	}
	if (lowHz < 4000) {
		// 高域の情報がないので解析しない
		return;
	}

	swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 200;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
	highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	if (swidth == 0) {
		swidth = 1;
	}

	// パワーの平均を計算する
	//
	avg = 0;
	for (i = lowIndex,n = 0;i < highIndex;i++) {
		avg += ovInfo->power[i];
		n++;
	}
	if (n > 0) {
		avg /= n;
	}
	avgLine = avg;
	ofs = lowIndex;
	//
	// powerが強いものを採用する(sig1)
	if (param->sig1Enb == 1) {
		window = swidth;
		maxOfs = ofs;
		maxWin = swidth;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			maxAvg = -1;
			skipCnt = 0;
			for (window = swidth;window < width;window++) {
				skipCnt++;
				if (param->hfaFast && (skipCnt % 4) != 0) {
					continue;
				}
				avg = 0;
				for (i = ofs,n = 0;i < highIndex;i += window) {
					avg += ovInfo->power[i];
					n++;
				}
				if (n > 0) {
					avg /= n;
				}

				if (maxAvg == -1 || maxAvg < avg) {
					maxAvg = avg;
					maxOfs = ofs;
					maxWin = window;
				}

			}
			if (avgLine * param->sig1AvgLineNx < maxAvg) {
				pha = ovInfo->phase[maxOfs];
				for (i = maxOfs,n = 1;i < nSample;i += maxWin,n++) {
					if (ovInfo->pw[i] < maxAvg / n) {
						ovInfo->pw[i] = maxAvg / n;
						ovInfo->phase[i] = pha;

						pha += param->sig1Phase;
						if (pha >= 360) {
							pha -= 360;
						}
						if (pha < 0) {
							pha += 360;
						}
					}
				}
			}
		}
	}

	//
	// window間隔の前後のpowerの差の累積が少ないものを採用する(sig2)
	if (param->sig2Enb == 1) {
		for (ofs = lowIndex;ofs < lowIndex + width;ofs++) {
			maxAvg = -1;
			skipCnt = 0;
			for (window = swidth;window < width;window++) {
				skipCnt++;
				if (param->hfaFast && (skipCnt % 4) != 0) {
					continue;
				}
				avg = 0;
				for (i = ofs,n = 0;i < highIndex;i += window) {
					if (1 + window < i) {
						if (ovInfo->power[i - window] <= ovInfo->power[i]) {
							avg += ovInfo->power[i] - ovInfo->power[i - window];
						} else {
							avg += ovInfo->power[i] - ovInfo->power[i - window];
						}
						n++;
					}
				}
				if (n > 0) {
					avg /= n;
				}

				if (maxAvg == -1 || maxAvg > avg) {
					maxAvg = avg;
					maxOfs = ofs;
					maxWin = window;
				}

			}
			avg = 0;
			for (i = maxOfs,n = 0;i < highIndex;i += maxWin) {
				avg += ovInfo->power[i];
				n++;
			}
			if (n > 0) {
				avg /= n;
			}
			maxAvg = avg;
			pha = ovInfo->phase[maxOfs];
			for (i = maxOfs,n = 1;i < nSample;i += maxWin,n++) {
				if (ovInfo->pw[i] > maxAvg) {
					ovInfo->pw[i] = maxAvg;
					ovInfo->phase[i] = pha;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha -= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
				} else if (ovInfo->pw[i] == 0) {
					ovInfo->pw[i] = maxAvg / 10;
					ovInfo->phase[i] = pha;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha -= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
				}
			}
		}
	}
}
#if 1	// 0.7 の HFA3
void anaOverToneHFA3(OVERTONE_INFO *ovInfo,PARAM_INFO *param)
{
	DWORD ofs,window,width,swidth;
	DWORD validIndex;
	DWORD n,i,j,k;
	double avg;
	double avgRef;
	double avgLine;
	double diff,diff0,diff1,diff2,diff3,diff4,diff5,diffP;
	double refPw[8];
	double avgPw,avgPw2,avgPw3;
	double tmpAvgPw,tmpAvgPw2,tmpAvgPw3;
	int    avgPwNX,avgPwNX2,avgPwNX3;
	long   step;
	long skipCnt;
	double tbl_hfaDiffMin[5] = {0.84,1.04,1.24,1.74,2.14};

	// 予測との最小パラメータ保存
	int minWin;
	int minType;
	int max_i;
	double minDiff;
	int nn;
	int odd;
	double hz;
	DWORD  baseOfs;
	double tmpPw,tmpPw2;
	int lowIndex,highIndex;
	int lowRange,highRange;
	int minWidth;
	int pha;
	double phaRand;
	int phaTmp = 0;
	long nSample;
	long lowHz,wid;
	double areaAvg;
	nSample = ovInfo->nSample;

	//
	// 初期化
	for (i = 1;i < ovInfo->nSample;i++) {
		ovInfo->pw[i] = 0;
		ovInfo->diff[i] = -1;
		ovInfo->base[i] = 0;
		ovInfo->baseToAdj[i] = 0;
		ovInfo->sign[i] = 0;
	}
	swidth = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1600;
#if 1
	for (i = 1;i < ovInfo->nSample;i+= swidth) {
		avg = 0;
		for (j = i,n = 0;n < swidth && j < ovInfo->nSample;n++) {
			avg += ovInfo->power[j];
		}
		avg /= n;
		for (j = i,n = 0;n < swidth && j < ovInfo->nSample;n++) {
			ovInfo->baseToAdj[j] = ovInfo->power[j] - avg;
			ovInfo->power[j] -= ovInfo->baseToAdj[j];
			ovInfo->base[j] = avg;
		}
	}
#endif
	if (ovInfo->validSamplingRate < 8000) {
		// 高域の情報がないので解析しない
		return;
	}

	lowHz	= 8000;
	wid		= 3000;
	if (lowHz + wid  >= ovInfo->validSamplingRate) {
		wid = ovInfo->validSamplingRate - lowHz;
	}

	step	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 50;
	step	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 77;	// 20181007
	step	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 97;	// 20191024
	step	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 120;	// 20191024
	step	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 248;	// 20191024
	swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 190;
	swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 177;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
	minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
	if (ovInfo->validSamplingRate < 16000) {
		highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	} else {
		highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
	}
	if (swidth == 0) {
		swidth = 1;
	}
	if (step == 0) {
		step = 1;
	}
	avg = 0;
	for (i = lowIndex,n = 0;i < highIndex;i++) {
		avg += ovInfo->power[i];
		n++;
	}
	if (n > 0) {
		avg /= n;
	}
	avgLine = avg;

	if (param->sig2Enb == 1) {
		// 前後のwindowで振幅の差が少ない音声の補間
		window = width;
		minWin = window;
		minType = 0;

		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 50;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 77;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 97;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 120;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 248;
			for (window = swidth;window < width;window+=step) {
//				step = 1;
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = diff1 = diff2 = diff3 = diff4 = diff5 = diffP = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				odd = 1;	// 奇数倍のみ倍音がある
				refPw[0] = -1;
				refPw[1] = -1;
				refPw[2] = -1;
				refPw[3] = -1;
				refPw[4] = -1;
				refPw[5] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				ovInfo->sign[baseOfs] = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						refPw[1] = ovInfo->power[i] * n;
						refPw[2] = ovInfo->power[i] * odd;
						refPw[3] = ovInfo->power[i] * (odd * odd);
						refPw[4] = ovInfo->power[i];
						refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					}
					// 平均より大きい音声か、小さい音声か
					if (ovInfo->power[i] > ovInfo->base[i]) {
						ovInfo->sign[baseOfs]++;
					} else if (ovInfo->power[i] < ovInfo->base[i]) {
						ovInfo->sign[baseOfs]--;
					}
					
					// 前後のパワーの差の計算
					if (i - window >= ofs) {
						if (ovInfo->power[i - window] >= ovInfo->power[i]) {
							diff = ovInfo->power[i - window] - ovInfo->power[i];
						} else {
							diff = ovInfo->power[i] - ovInfo->power[i - window];
						}
					}
					diffP += (diff * tbl_hfaDiffMin[param->hfaDiffMin - 1]);

					avgPw += ovInfo->power[i];
					avgPwNX++;
					if ((avgPwNX & 0x01) == 0) {
						avgPw2 += ovInfo->power[i];
						avgPwNX2++;
					}
					if ((avgPwNX % 3) == 0) {
						avgPw3 += ovInfo->power[i];
						avgPwNX3++;
					}
					
					// 1/f(振幅が1/fになっているもの)
					diff = refPw[0] / hz;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff0 += diff;

					// 鋸波(nの逆数で小さくなる)
					diff = refPw[1] / n;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff1 += diff;

					// 短形波(奇数倍音,nの逆数で小さくなる)
					diff = refPw[2] / odd;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff2 += diff;

					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					diff = refPw[3] / (odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff3 += diff;

					// パルス(n番目の倍音でもパワーは同じ)
					diff = refPw[4];
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff4 += diff;

					// その他(もっとパワーが小さくなるパターン)
					diff = refPw[5] / (odd * odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff5 += diff;

					nn++;
				}

//				diff0 += diffP * 2;
//				diff1 += diffP * 2;
//				diff2 += diffP * 2;
//				diff3 += diffP * 2;
//				diff4 += diffP * 2;
//				diff5 += diffP * 2;

				if (nn > 0) {
					diff0 /= nn;
					diff1 /= nn;
					diff2 /= nn;
					diff3 /= nn;
					diff4 /= nn;
					diff5 /= nn;
					diffP /= nn;
//					if (refPw[4] > 0) {
//						diff0 /= refPw[4];
//						diff1 /= refPw[4];
//						diff2 /= refPw[4];
//						diff3 /= refPw[4];
//						diff4 /= refPw[4];
//						diff5 /= refPw[4];
//						diffP /= refPw[4];
//					}
				}

				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
					tmpAvgPw  = avgPw / avgPwNX;
					tmpAvgPw2 = avgPw2 / avgPwNX2;
					tmpAvgPw3 = avgPw3 / avgPwNX3;
					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {						continue;
					}
				}

				if (minDiff == -1 || minDiff > diffP) {
					minDiff = diffP;
					minWin = window;
					minType = 0;
				}
				if (minDiff > diff1) {
					minDiff = diff1;
					minWin = window;
					minType = 1;
				}
				if (minDiff > diff2) {
					minDiff = diff2;
					minWin = window;
					minType = 2;
				}
				if (minDiff > diff3) {
					minDiff = diff3;
					minWin = window;
					minType = 3;
				}
				if (minDiff > diff4) {
					minDiff = diff4;
					minWin = window;
					minType = 4;
				}
				if (minDiff > diff5) {
					minDiff = diff5;
					minWin = window;
					minType = 5;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			refPw[4] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}

			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					refPw[1] = ovInfo->power[i] * n;
					refPw[2] = ovInfo->power[i] * odd;
					refPw[3] = ovInfo->power[i] * (odd * odd);
					refPw[4] = ovInfo->power[i];
					refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (minType == 0) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[0] / hz;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.41;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.41;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->sign[baseOfs] > 8) {
							if (tmpPw2 < tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						} else if (ovInfo->sign[baseOfs] < 0) {
							if (tmpPw2 > tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						}
					}
				} else if (minType == 1) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[1] / n;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.41;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.41;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->sign[baseOfs] > 8) {
							if (tmpPw2 < tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						} else if (ovInfo->sign[baseOfs] < 0) {
							if (tmpPw2 > tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						}
					}
				} else if (minType == 2) {
					// 短形波(奇数倍音,nの逆数で小さくなる)
					tmpPw = refPw[2] / odd;
					phaRand = 1;
					pha = ovInfo->phase[max_i];
					phaTmp = pha;
					if (n & 0x01) {
						phaTmp = pha + 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.41;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.41;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->sign[baseOfs] > 8) {
							if (tmpPw2 < tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						} else if (ovInfo->sign[baseOfs] < 0) {
							if (tmpPw2 > tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						}
					}
				} else if (minType == 3) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					phaTmp = pha;
					if (n & 0x01) {
						phaTmp = pha + 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.41;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.41;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->sign[baseOfs] > 8) {
							if (tmpPw2 < tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						} else if (ovInfo->sign[baseOfs] < 0) {
							if (tmpPw2 > tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						}
					}
				} else if (minType == 4) {
					// パルス(n番目の倍音でもパワーは同じ)
					tmpPw = refPw[4];
					phaRand = rand() * 6;
					phaRand -= 3;
					pha += phaRand;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.18;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.18;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->sign[baseOfs] > 8) {
							if (tmpPw2 < tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.18;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						} else if (ovInfo->sign[baseOfs] < 0) {
							if (tmpPw2 > tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.18;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						}
					}
				} else if (minType == 5) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.41;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.41;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->sign[baseOfs] > 8) {
							if (tmpPw2 < tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						} else if (ovInfo->sign[baseOfs] < 0) {
							if (tmpPw2 > tmpPw) {
								ovInfo->pw[i] = tmpPw * 0.41;
								ovInfo->phase[i] = pha;
								ovInfo->diff[i] = minDiff;
							}
						}
					}
				}
			}
		}
	}

#if 1
	if (param->sig1Enb == 1) {
		// powerが強いもの優先で補間する
		lowHz	= 9500;
		wid		= 2500;
		if (lowHz + wid  >= ovInfo->validSamplingRate) {
			wid = ovInfo->validSamplingRate - lowHz;
		}

		swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 190;
		width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
		lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
		minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
		if (ovInfo->validSamplingRate < 16000) {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
		} else {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
		}
		if (swidth == 0) {
			swidth = 1;
		}

		areaAvg = 0;
		for (i = lowIndex,n = 0;n < width;i++,n++) {
			areaAvg += ovInfo->power[i];
		}
		if (n > 0) {
			areaAvg /= n;
		}

		window = width;
		minWin = window;
		minType = 0;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			if (ovInfo->power[ofs] < areaAvg) {
				continue;
			}
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 50;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 103;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 77;
			for (window = swidth;window < (width / 2);window+=step) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}
				//step = 1;
				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				ovInfo->sign[baseOfs] = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					// 平均より大きい音声か、小さい音声か
					if (ovInfo->power[i] > ovInfo->base[i]) {
						ovInfo->sign[baseOfs]++;
					} else if (ovInfo->power[i] < ovInfo->base[i]) {
						ovInfo->sign[baseOfs]--;
					}

					if (i - window >= ofs) {
						if (ovInfo->power[i - window] < ovInfo->power[i]) {
							avgPw /= 75;
						}
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						max_i = i;
					}
					if (ovInfo->power[i] > areaAvg) {
						avgPw += ovInfo->power[i];
					} else {
						avgPw /= 75;
					}
					avgPwNX++;
					nn++;
				}

				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
					tmpAvgPw  = avgPw / avgPwNX;
					tmpAvgPw2 = avgPw2 / avgPwNX2;
					tmpAvgPw3 = avgPw3 / avgPwNX3;
					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {						continue;
					}
				}
				if (avgPwNX > 0) {
					avgPw /= avgPwNX;
				}
				if (minDiff == -1 || minDiff < avgPw) {
					minDiff = avgPw;
					minWin = window;
					minType = 0;
				}
			}

			// 一番累積のパワーが強いものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				// 1/f(振幅が1/fになっているもの)
				tmpPw = refPw[0] / hz;
				phaRand = 1;
				pha += param->sig1Phase;
				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
				}
				if (ovInfo->pw[i] == 0 && ovInfo->sign[baseOfs] > 16) {
					ovInfo->pw[i] = tmpPw * 0.55;
					ovInfo->phase[i] = pha;
					ovInfo->diff[i] = -1;
				}
			}
		}
	}
#endif

	if (param->sig3Enb == 1) {
//	if (0) {
		// powerが弱いものを補間する
		lowHz	= 8000;
		wid		= 4000;
		if (lowHz + wid  >= ovInfo->validSamplingRate) {
			wid = ovInfo->validSamplingRate - lowHz;
		}

		swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 180;
		width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
		lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
		minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
		if (ovInfo->validSamplingRate < 16000) {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
		} else {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
		}
		if (swidth == 0) {
			swidth = 1;
		}

		window = width;
		minWin = swidth;
		
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 50;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 66;
			step = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 77;
			for (window = swidth;window < (width / 2);window+=step) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}
				step = 1;
				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				ovInfo->sign[baseOfs] = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}
					// 平均より大きい音声か、小さい音声か
					if (ovInfo->power[i] > ovInfo->base[i]) {
						ovInfo->sign[baseOfs]++;
					} else if (ovInfo->power[i] < ovInfo->base[i]) {
						ovInfo->sign[baseOfs]--;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						max_i = i;
					}
					avgPw += ovInfo->power[i];
					avgPwNX++;
					nn++;
				}

				if (avgPwNX > 0) {
					avgPw /= avgPwNX;
				}
				if (minDiff == -1 || minDiff > avgPw) {
					minDiff = avgPw;
					minWin = window;
					minType = 0;
				}
			}

			// 一番累積のパワーが弱いものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				// 1/f(振幅が1/fになっているもの)
				tmpPw = refPw[0] / hz;
				phaRand = 1;
				pha += param->sig3Phase;
				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
				}
				if (ovInfo->pw[i] == 0 && ovInfo->sign[baseOfs] < -16) {
					ovInfo->pw[i] = tmpPw * 0.18;
					ovInfo->phase[i] = pha;
					ovInfo->diff[i] = -1;
				}
			}
		}
	}

	if (param->sig2Enb == 1 && param->hfaWide) {
		lowHz	= 7000;
		wid		= 4000;
		if (lowHz + wid  >= ovInfo->validSamplingRate) {
			wid = ovInfo->validSamplingRate - lowHz;
		}

		swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 180;
		width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
		lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
		minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
		if (ovInfo->validSamplingRate < 16000) {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
		} else {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
		}
		if (swidth == 0) {
			swidth = 1;
		}

		// 前後のwindowで振幅の差が少ない音声の補間
		window = width;
		minWin = window;
		minType = 0;

		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window+=step) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = diff1 = diff2 = diff3 = diff4 = diff5 = diffP = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				odd = 1;	// 奇数倍のみ倍音がある
				refPw[0] = -1;
				refPw[1] = -1;
				refPw[2] = -1;
				refPw[3] = -1;
				refPw[4] = -1;
				refPw[5] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						refPw[1] = ovInfo->power[i] * n;
						refPw[2] = ovInfo->power[i] * odd;
						refPw[3] = ovInfo->power[i] * (odd * odd);
						refPw[4] = ovInfo->power[i];
						refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					}
					// 前後のパワーの差の計算
					if (i - window >= ofs) {
						if (ovInfo->power[i - window] >= ovInfo->power[i]) {
							diff = ovInfo->power[i - window] - ovInfo->power[i];
						} else {
							diff = ovInfo->power[i] - ovInfo->power[i - window];
						}
					}
					diffP += (diff * 1.48);

					avgPw += ovInfo->power[i];
					avgPwNX++;
					if ((avgPwNX & 0x01) == 0) {
						avgPw2 += ovInfo->power[i];
						avgPwNX2++;
					}
					if ((avgPwNX % 3) == 0) {
						avgPw3 += ovInfo->power[i];
						avgPwNX3++;
					}
					
					// 1/f(振幅が1/fになっているもの)
					diff = refPw[0] / hz;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff0 += diff;

					// 鋸波(nの逆数で小さくなる)
					diff = refPw[1] / n;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff1 += diff;

					// 短形波(奇数倍音,nの逆数で小さくなる)
					diff = refPw[2] / odd;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff2 += diff;

					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					diff = refPw[3] / (odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff3 += diff;

					// パルス(n番目の倍音でもパワーは同じ)
					diff = refPw[4];
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff4 += diff;

					// その他(もっとパワーが小さくなるパターン)
					diff = refPw[5] / (odd * odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff5 += diff;

					nn++;
				}

				diff0 += diffP * 4;
				diff1 += diffP * 4;
				diff2 += diffP * 4;
				diff3 += diffP * 4;
				diff4 += diffP * 4;
				diff5 += diffP * 4;

				if (nn > 0) {
					diff0 /= nn;
					diff1 /= nn;
					diff2 /= nn;
					diff3 /= nn;
					diff4 /= nn;
					diff5 /= nn;
					diffP /= nn;
				}

				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
					tmpAvgPw  = avgPw / avgPwNX;
					tmpAvgPw2 = avgPw2 / avgPwNX2;
					tmpAvgPw3 = avgPw3 / avgPwNX3;
					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {						continue;
					}
				}

				if (minDiff == -1 || minDiff > diffP) {
					minDiff = diffP;
					minWin = window;
					minType = 0;
				}
				if (minDiff > diff1) {
					minDiff = diff1;
					minWin = window;
					minType = 1;
				}
				if (minDiff > diff2) {
					minDiff = diff2;
					minWin = window;
					minType = 2;
				}
				if (minDiff > diff3) {
					minDiff = diff3;
					minWin = window;
					minType = 3;
				}
				if (minDiff > diff4) {
					minDiff = diff4;
					minWin = window;
					minType = 4;
				}
				if (minDiff > diff5) {
					minDiff = diff5;
					minWin = window;
					minType = 5;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			refPw[4] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}

			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					refPw[1] = ovInfo->power[i] * n;
					refPw[2] = ovInfo->power[i] * odd;
					refPw[3] = ovInfo->power[i] * (odd * odd);
					refPw[4] = ovInfo->power[i];
					refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (minType == 0) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[0] / hz;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.77 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.32;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 1) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[1] / n;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.77 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.28;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 2) {
					// 短形波(奇数倍音,nの逆数で小さくなる)
					tmpPw = refPw[2] / odd;
					phaRand = 1;
					pha = ovInfo->phase[max_i];
					phaTmp = pha;
					if (n & 0x01) {
						phaTmp = pha + 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.77 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.28;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 3) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					phaTmp = pha;
					if (n & 0x01) {
						phaTmp = pha + 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.77 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.28;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 4) {
					// パルス(n番目の倍音でもパワーは同じ)
					tmpPw = refPw[4];
					phaRand = rand() * 6;
					phaRand -= 3;
					pha += phaRand;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.77 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.28;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 5) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.77 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.28;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					}
				}
			}
		}
	}

//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			for (i = 300;i < 750;i++) {
//				fprintf(ofp,"%lf,",ovInfo->pw[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//	}

//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			for (i = 300;i < 750;i++) {
//				fprintf(ofp,"%lf,",ovInfo->pw[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//	}

#if 0
	//
	// 位相の調整
	for (i = 1;i < validIndex;i++) {
		ovInfo->diff[i] = -1;
	}
	window = width;
	for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
		minDiff = -1;
		skipCnt = 0;
		for (window = swidth;window < width;window++) {
			if (window < minWidth) {
				if ((ofs - lowIndex) > window * 1) {
					continue;
				}
			} else {
				if ((ofs - lowIndex) > window) {
					continue;
				}
			}
			skipCnt++;
			if (param->hfaFast && (skipCnt % 8) != 0) {
				continue;
			}
			// 位相を調べる
			diffP = 0;
			baseOfs = ofs - ((ofs / window) * window);
			if (baseOfs == 0) {
				baseOfs = window;
			}
			n = 1;
			refPw[0] = -1;
			for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}

				if (refPw[0] == -1) {
					refPw[0] = ovInfo->phase[i];
				}
				// 前後の位相の差の計算
				if (i - window >= ofs) {
					if (refPw[0] <= ovInfo->phase[i]) {
						diff = ovInfo->phase[i] - refPw[0];
					} else {
						diff = refPw[0] - ovInfo->phase[i];
					}
				}
				diffP += diff;
				nn++;
			}

			if (nn > 0) {
				diffP /= nn;
			}

			if (minDiff == -1 || minDiff > diffP) {
				minDiff = diffP;
				minWin = window;
			}
		}

		// 一番予測誤差が少なかったものを採用する。

		baseOfs = ofs - ((ofs / minWin) * minWin);
		if (baseOfs == 0) {
			baseOfs = minWin;
		}

		pha = ovInfo->phase[baseOfs];
		n = 1;		// 奇数偶数倍すべてに倍音がある

		refPw[0] = -1;

		for (i = baseOfs;i < validIndex;i += minWin,n++) {
			hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
			if (hz < lowHz) {
				continue;
			}
			if (refPw[0] == -1) {
				refPw[0] = ovInfo->phase[i];
			}

			phaRand = rand() * 2;
			phaRand -= 1;
			pha = refPw[0];
			//pha += phaRand;
			if (pha >= 360) {
				pha %= 360;
			}
			if (pha < 0) {
				pha += 360;
			}
			if (ovInfo->diff[i] == -1 || ovInfo->diff[i] > minDiff) {
				ovInfo->phase[i] = pha;
				ovInfo->diff[i] = minDiff;
			}
		}
	}
#endif
	// 補間されていない箇所の対応
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	if (param->sig3Phase >= -6 || param->sig3Phase <= 6) {
		// 位相の修正
		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->phase[i];
					}
					
					diff = refPw[0];
					if (diff >= ovInfo->phase[i]) {
						diff = diff - ovInfo->phase[i];
					} else {
						diff = ovInfo->phase[i] - diff;
					}
					diff0 += diff;
					nn++;
				}

				if (nn > 0) {
					diff0 /= nn;
				}

				if (minDiff == -1 || minDiff > diff0) {
					minDiff = diff0;
					minWin = window;
					minType = 0;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある

			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
					pha += 360;
					pha += 360;
					pha %= 360;
				}
				if (ovInfo->pw[i] != 0) {
					if (ovInfo->diff[i] == -1 || ovInfo->diff[i] > minDiff) {
						ovInfo->phase[i] = pha;
					}
				}
			}
		}
	}
#if 1
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * (ovInfo->validSamplingRate - 5000);
	highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * (ovInfo->validSamplingRate - 1000);
	i = 1;
	do {
		for (j = i,n = 0,k = lowIndex;n < highIndex - lowIndex && j < ovInfo->nSample;j++,n++) {
			ovInfo->pw[j] += ovInfo->baseToAdj[k];
			k++;
			if (k > highIndex) {
				k = lowIndex;
			}
		}
		i += n;
	} while (i < ovInfo->nSample);
#endif
#if 0
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	highIndex = nSample;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 5000;
	i = lowIndex;
	do {
		avg = 0;
		for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
			avg += ovInfo->pw[j];
		}
		if (n > 0) {
			avg /= n;
			avg = log10(avg) * 20;
		}
		for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
			tmpPw = ovInfo->pw[j];
			if (tmpPw > 0) {
				tmpPw = log10(tmpPw) * 20;
			}
			if (tmpPw + 15 < avg) {
				ovInfo->pw[j] *= 13;
			} else if (tmpPw > avg + 15) {
				ovInfo->pw[j] /= 13;
			}
		}
		i = j;
	} while (i < ovInfo->nSample);
#endif
//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			for (i = 0;i < 4096;i++) {
//				fprintf(ofp,"%lf,%lf\n",ovInfo->pw[i],ovInfo->phase[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//		ovInfo->log = 0;
//	}

#if 1
	if (param->hfaWide) {
		// 解析情報を平均化し、急激な変化を抑える
		for (i = baseOfs;i + 1< nSample;i++) {
			hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
			if (hz < lowHz) {
				continue;
			}
			if (ovInfo->pw[i] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
#endif
	if (1) {
		double nx;
		// 付加した広域で強すぎるものがあれば弱める。
		lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
		highIndex = nSample;
		swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 5000;
		i = lowIndex;
		nx = 6.7;
		if (ovInfo->validSamplingRate <= 12000) {
			nx = 3.2;
		}
		do {
			avg = 0;
			for (j = i,n = 0;n < swidth && j < ovInfo->nSample;j++,n++) {
				avg += ovInfo->pw[j];
			}
			if (n > 0) {
				avg /= n;
			}
			for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
				if (avg * nx < ovInfo->pw[j]) {
					ovInfo->pw[j] *= avg * nx / ovInfo->pw[j];
				}
			}
			i = j;
		} while (i < ovInfo->nSample);
	}
if (0) {
	if (ovInfo->pw[i] > 0 && ovInfo->pw[i + 1] > 0) {
		ovInfo->pw[i] = (ovInfo->pw[i] + ovInfo->pw[i+1]) / 2;
		ovInfo->phase[i] = (ovInfo->phase[i] + ovInfo->phase[i+1]) / 2;
		if (ovInfo->phase[i] >= 360) {
			pha = ovInfo->phase[i];
			ovInfo->phase[i] = (pha % 360);
		}
		if (ovInfo->phase[i] < 0) {
			ovInfo->phase[i] += 360;
		}
	}
}
}
#else
void anaOverToneHFA3(OVERTONE_INFO *ovInfo,PARAM_INFO *param)
{
	DWORD ofs,window,width,swidth;
	DWORD n,i,j,k,lll;
	int ii,jj;
	DWORD validN;
	double avg,minAvg,maxAvg,peekPw,diffPw,diffAvg;
	double keyPw;
	double avgLine;
	double diff,diff0,diff1,diff2,diff3,diff4,diff5,diff6,diff7,diffP;
	double refPw[8];
	double avgPw,avgPw2,avgPw3;
	double tmpAvgPw,tmpAvgPw2,tmpAvgPw3;
	int    avgPwNX,avgPwNX2,avgPwNX3;
	double maxPwAvg;
	long skipCnt;

	// 予測との最小パラメータ保存
	int minWin;
	int minType;
	int max_i;
	double minDiff;

	int nd,nn;
	int odd;
	double ndLv;
	double hz,hz2;
	DWORD  baseOfs;
	int minOfs,maxOfs,maxWin,maxType;
	double tmpPw;
	double nx,nx2;
	int tmpOfs,tmpWin;
	int lowIndex,highIndex;
	int minWidth;
	int miOffset;
	int pha;
	double phaRand;
	int phaTmp = 0;
	long nSample;
	long lowHz,wid;
	double areaAvg;
	nSample = ovInfo->nSample;

	//
	// 初期化
	for (i = 1;i < ovInfo->nSample;i++) {
		ovInfo->pw[i] = 0;
		ovInfo->diff[i] = -1;
		ovInfo->baseToAdj[i] = 0;
	}
	swidth = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 2000;
	//swidth = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 400;
#if 0
	for (i = 1;i < ovInfo->nSample;i+= swidth) {
		avg = 0;
		for (j = i,n = 0;n < swidth && j < ovInfo->nSample;n++) {
			avg += ovInfo->power[j];
		}
		avg /= n;
		for (j = i,n = 0;n < swidth && j < ovInfo->nSample;n++) {
			ovInfo->baseToAdj[j] = ovInfo->power[j] - avg;
			ovInfo->power[j] -= ovInfo->baseToAdj[j];
		}
	}
#endif
	lowHz	= 6000;
	wid		= 3500;
	if (param->hfaWide) {
		lowHz	= 4500;
		wid		= 5500;
	}
	if (lowHz + wid + 2000 >= ovInfo->validSamplingRate) {
		wid = 2000;
		lowHz = 4500;
	}

	if (lowHz < 4000) {
		// 高域の情報がないので解析しない
		return;
	}

	swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 180;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
	minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
	if (ovInfo->validSamplingRate < 16000) {
		highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	} else {
		highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
	}
	if (swidth == 0) {
		swidth = 1;
	}

	avg = 0;
	for (i = lowIndex,n = 0;i < highIndex;i++) {
		avg += ovInfo->power[i];
		n++;
	}
	if (n > 0) {
		avg /= n;
	}
	avgLine = avg;

	if (param->sig2Enb == 1) {
		// 前後のwindowで振幅の差が少ない音声の補間
		window = width;
		minWin = window;
		minType = 0;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = diff1 = diff2 = diff3 = diff4 = diff5 = diffP = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				odd = 1;	// 奇数倍のみ倍音がある
				refPw[0] = -1;
				refPw[1] = -1;
				refPw[2] = -1;
				refPw[3] = -1;
				refPw[4] = -1;
				refPw[5] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						refPw[1] = ovInfo->power[i] * n;
						refPw[2] = ovInfo->power[i] * odd;
						refPw[3] = ovInfo->power[i] * (odd * odd);
						refPw[4] = ovInfo->power[i];
						refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					}
					// 前後のパワーの差の計算
					if (i - window >= ofs) {
						if (ovInfo->power[i - window] >= ovInfo->power[i]) {
							diff = ovInfo->power[i - window] - ovInfo->power[i];
						} else {
							diff = ovInfo->power[i] - ovInfo->power[i - window];
						}
					}
					diffP += diff;

					avgPw += ovInfo->power[i];
					avgPwNX++;
					if ((avgPwNX & 0x01) == 0) {
						avgPw2 += ovInfo->power[i];
						avgPwNX2++;
					}
					if ((avgPwNX % 3) == 0) {
						avgPw3 += ovInfo->power[i];
						avgPwNX3++;
					}
					
					// 1/f(振幅が1/fになっているもの)
					diff = refPw[0] / hz;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff0 += diff;

					// 鋸波(nの逆数で小さくなる)
					diff = refPw[1] / n;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff1 += diff;

					// 短形波(奇数倍音,nの逆数で小さくなる)
					diff = refPw[2] / odd;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff2 += diff;

					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					diff = refPw[3] / (odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff3 += diff;

					// パルス(n番目の倍音でもパワーは同じ)
					diff = refPw[4];
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff4 += diff;

					// その他(もっとパワーが小さくなるパターン)
					diff = refPw[5] / (odd * odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff5 += diff;

					nn++;
				}

				if (nn > 0) {
					diff0 /= nn;
					diff1 /= nn;
					diff2 /= nn;
					diff3 /= nn;
					diff4 /= nn;
					diff5 /= nn;
					diffP /= nn;
				}

				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
					tmpAvgPw  = avgPw / avgPwNX;
					tmpAvgPw2 = avgPw2 / avgPwNX2;
					tmpAvgPw3 = avgPw3 / avgPwNX3;
					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {						continue;
					}
				}

				if (minDiff == -1 || minDiff > diffP) {
					minDiff = diffP;
					minWin = window;
					minType = 0;
				}
				if (minDiff > diff1) {
					minDiff = diff1;
					minWin = window;
					minType = 1;
				}
				if (minDiff > diff2) {
					minDiff = diff2;
					minWin = window;
					minType = 2;
				}
				if (minDiff > diff3) {
					minDiff = diff3;
					minWin = window;
					minType = 3;
				}
				if (minDiff > diff4) {
					minDiff = diff4;
					minWin = window;
					minType = 4;
				}
				if (minDiff > diff5) {
					minDiff = diff5;
					minWin = window;
					minType = 5;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			refPw[4] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					refPw[1] = ovInfo->power[i] * n;
					refPw[2] = ovInfo->power[i] * odd;
					refPw[3] = ovInfo->power[i] * (odd * odd);
					refPw[4] = ovInfo->power[i];
					refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (minType == 0) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[0] / hz;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 1) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[1] / n;
					phaRand = 1;
					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 2) {
					// 短形波(奇数倍音,nの逆数で小さくなる)
					tmpPw = refPw[2] / odd;
					phaRand = 1;
					pha = ovInfo->phase[max_i];
					phaTmp = pha + param->sig2Phase;
					if (n & 0x01) {
						phaTmp += 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 3) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					phaTmp = pha + param->sig2Phase;
					if (n & 0x01) {
						phaTmp += 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 4) {
					// パルス(n番目の倍音でもパワーは同じ)
					tmpPw = refPw[4];
					phaRand = rand() * 6;
					phaRand -= 3;
					pha += phaRand;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 5) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				}
			}
		}
	}
//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			for (i = 300;i < 750;i++) {
//				fprintf(ofp,"%lf,",ovInfo->pw[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//	}

#if 0
	if (param->sig3Enb == 1 && (param->sig3Phase < -6 || param->sig3Phase > 6)) {
		// 位相の修正
		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->phase[i];
					}
					
					diff = refPw[0];
					if (diff >= ovInfo->phase[i]) {
						diff = diff - ovInfo->phase[i];
					} else {
						diff = ovInfo->phase[i] - diff;
					}
					diff0 += diff;
					nn++;
				}

				if (nn > 0) {
					diff0 /= nn;
				}

				if (minDiff == -1 || minDiff > diff0) {
					minDiff = diff0;
					minWin = window;
					minType = 0;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある

			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
					pha += 360;
					pha += 360;
					pha %= 360;
				}
				if (ovInfo->pw[i] != 0) {
					if (ovInfo->diff[i] == -1 || ovInfo->diff[i] > minDiff) {
						ovInfo->phase[i] = pha;
					}
				}
			}
		}
	}
#endif
	if (param->sig1Enb == 1) {
		// powerが強いもの優先で補間する

		areaAvg = 0;
		for (i = lowIndex,n = 0;n < width;i++,n++) {
			areaAvg += ovInfo->power[i];
		}
		if (n > 0) {
			areaAvg /= n;
		}

		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			if (ovInfo->power[ofs] < areaAvg) {
				continue;
			}
			
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}
				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}


					if (i - window >= ofs) {
						if (ovInfo->power[i - window] < ovInfo->power[i]) {
							avgPw /= 75;
						}
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						max_i = i;
					}
					if (ovInfo->power[i] > areaAvg) {
						avgPw += ovInfo->power[i];
					} else {
						avgPw /= 75;
					}
					avgPwNX++;
					nn++;
				}

//				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
//					tmpAvgPw  = avgPw / avgPwNX;
//					tmpAvgPw2 = avgPw2 / avgPwNX2;
//					tmpAvgPw3 = avgPw3 / avgPwNX3;
//					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {//						continue;
//					}
//				}
				if (avgPwNX > 0) {
					avgPw /= avgPwNX;
				}
				if (minDiff == -1 || minDiff < avgPw) {
					minDiff = avgPw;
					minWin = window;
					minType = 0;
				}
			}

			// 一番累積のパワーが強いものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある

			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				// 1/f(振幅が1/fになっているもの)
				tmpPw = refPw[0] / hz;
				phaRand = 1;
				pha += param->sig1Phase;
				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
				}
				if (ovInfo->pw[i] == 0) {
					ovInfo->pw[i] = (tmpPw * 0.70);
					ovInfo->phase[i] = pha;
					ovInfo->diff[i] = -1;
				}
			}
		}
	}
//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			for (i = 300;i < 750;i++) {
//				fprintf(ofp,"%lf,",ovInfo->pw[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//	}

	if (param->sig3Enb == 1) {
		// powerが弱いものと補間値がないものを補間する

		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}
				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						max_i = i;
					}
					avgPw += ovInfo->power[i];
					avgPwNX++;
					nn++;
				}

				if (avgPwNX > 0) {
					avgPw /= avgPwNX;
				}
				if (minDiff == -1 || minDiff > avgPw) {
					minDiff = avgPw;
					minWin = window;
					minType = 0;
				}
			}

			// 一番累積のパワーが弱いものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				// 1/f(振幅が1/fになっているもの)
				tmpPw = refPw[0] / hz;
				phaRand = 1;
				pha += param->sig3Phase;
				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
				}
				if (ovInfo->pw[i] == 0) {
					ovInfo->pw[i] = (tmpPw * 0.35);
					ovInfo->phase[i] = pha;
					ovInfo->diff[i] = -1;
				}
			}
		}
	}

	// 補間されていない箇所の対応
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	if (param->sig3Phase >= -6 || param->sig3Phase <= 6) {
		// 位相の修正
		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->phase[i];
					}
					
					diff = refPw[0];
					if (diff >= ovInfo->phase[i]) {
						diff = diff - ovInfo->phase[i];
					} else {
						diff = ovInfo->phase[i] - diff;
					}
					diff0 += diff;
					nn++;
				}

				if (nn > 0) {
					diff0 /= nn;
				}

				if (minDiff == -1 || minDiff > diff0) {
					minDiff = diff0;
					minWin = window;
					minType = 0;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある

			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
					pha += 360;
					pha += 360;
					pha %= 360;
				}
				if (ovInfo->pw[i] != 0) {
					if (ovInfo->diff[i] == -1 || ovInfo->diff[i] > minDiff) {
						ovInfo->phase[i] = pha;
					}
				}
			}
		}
	}
#if 0
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * (ovInfo->validSamplingRate - 5000);
	highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * (ovInfo->validSamplingRate - 2000);
	i = 1;
	do {
		for (j = i,n = 0,k = lowIndex;n < highIndex - lowIndex && j < ovInfo->nSample;j++,n++) {
			ovInfo->pw[j] += ovInfo->baseToAdj[k];
			k++;
			if (k > highIndex) {
				k = lowIndex;
			}
		}
		i += n;
	} while (i < ovInfo->nSample);
#endif
#if 0
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	highIndex = nSample;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 5000;
	i = lowIndex;
	do {
		avg = 0;
		for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
			avg += ovInfo->pw[j];
		}
		if (n > 0) {
			avg /= n;
			avg = log10(avg) * 20;
		}
		for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
			tmpPw = ovInfo->pw[j];
			if (tmpPw > 0) {
				tmpPw = log10(tmpPw) * 20;
			}
			if (tmpPw + 15 < avg) {
				ovInfo->pw[j] *= 13;
			} else if (tmpPw > avg + 15) {
				ovInfo->pw[j] /= 13;
			}
		}
		i = j;
	} while (i < ovInfo->nSample);
#endif

#if 1
	if (param->hfaWide) {
		// 解析情報を平均化し、急激な変化を抑える
		for (i = baseOfs;i + 1< nSample;i++) {
			hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
			if (hz < lowHz) {
				continue;
			}
			if (ovInfo->pw[i] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
#endif

//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			for (i = 0;i < 4096;i++) {
//				fprintf(ofp,"%lf,%lf\n",ovInfo->pw[i],ovInfo->phase[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//		ovInfo->log = 0;
//	}
}
#endif
//---------------------------------------------------------------------------
// Function   : noiseCut
// Description: ノイズカット処理
// ---
//	nfs		 	:ノイズカットオフ周波数(この周波数以上の領域のノイズをカットする)
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	fp_w		:出力ファイル
//	param		:変換パラメータ
//
void noiseCut(long nfs,SSIZE inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3;
	long hfc;
	long outSampleR;
	long wkMemSize;
	long fftSize,i,j,n,nn,h;
	long lowIndex,highIndex;
	long hz,idx;
	double percent,per;
	int ignore_flag;
	double p;
	static double refPw[10000];
	double *pw;
	double cutLV[5]={0.99,1.06,1.12,1.18,1.30};
	SSIZE *pIn,*pOut;
	SSIZE startInSample,nSample;
	fftw_complex *fftw_in,*fftw_out;
	fftw_plan fftw_p,fftw_ip;

   	outSampleR = param->outSampleR;
	if (param->hfc != -1) {
		hfc = param->hfc;
	} else {
		hfc = param->inSampleR / 2;
	}
	if (hfc > param->inSampleR / 2) {
		hfc = param->inSampleR / 2;
	}

	fio_rewind(fp_r);
	fio_rewind(fp_w);

	if ((outSampleR == 32000) || (outSampleR == 44100) || (outSampleR == 48000)) {
		fftSize = 4096 / 2;
		fftSize = outSampleR / 14;
		fftSize = outSampleR / 8;
		fftSize = outSampleR / 10;
	}
	if ((outSampleR == 44100 * 2) || (outSampleR == 48000 * 2)) {
		fftSize = 8192 / 2;
		fftSize = outSampleR / 14;
		fftSize = outSampleR / 8;
		fftSize = outSampleR / 10;
	}
	if ((outSampleR == 32000 * 6) || (outSampleR == 44100 * 4) || (outSampleR == 48000 * 4)) {
		fftSize = 16384 / 2;
		fftSize = outSampleR / 14;
		fftSize = outSampleR / 8;
		fftSize = outSampleR / 10;
	}
	if ((outSampleR == 32000 * 12) || (outSampleR == 44100 * 8) || (outSampleR == 48000 * 8)) {
		fftSize = 32768 / 2;
		fftSize = outSampleR / 14;
		fftSize = outSampleR / 8;
		fftSize = outSampleR / 10;
	}
	if ((outSampleR == 32000 * 24) || (outSampleR == 44100 * 16) || (outSampleR == 48000 * 16)) {
		fftSize = 65536 / 2;
		fftSize = outSampleR / 14;
		fftSize = outSampleR / 8;
		fftSize = outSampleR / 10;
	}
	if ((outSampleR == 32000 * 48) || (outSampleR == 44100 * 32) || (outSampleR == 48000 * 32) || (outSampleR == 44100 * 64) || (outSampleR == 48000 * 64)) {
		fftSize = 65536;
		fftSize = outSampleR / 14;
		fftSize = outSampleR / 8;
		fftSize = outSampleR / 10;
	}

	wkMemSize = (65536 * 32) * sizeof (SSIZE);

	mem1 = (SSIZE *)malloc(wkMemSize);
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem2 = (SSIZE *)malloc(wkMemSize);
	if (mem2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem3 = (SSIZE *)malloc(wkMemSize);
	if (mem3 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	pw = (double *)malloc(65536*2 * sizeof (double));
	if (pw == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536*2);
	if (fftw_in == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536*2);
	if (fftw_out == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p = fftw_plan_dft_1d(fftSize,fftw_in,fftw_out,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip = fftw_plan_dft_1d(fftSize,fftw_out,fftw_in,FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	per = -1;
	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + (fftSize + (fftSize / 2));startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
			percent = ((double)startInSample / inSample);
			percent *= 100;
			percent = (double)((int)percent);
			if (percent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)percent);
				fflush(stdout);
				if (chkAbort(param,percent,2) == 1) exit(0);
			}
			per = percent;
		}

		memset(mem1,0,wkMemSize);

		if (startInSample >= 0 && startInSample + (fftSize * 12) < inSample) {
			nSample = fftSize * 12;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 12;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 12;
		}

		memset(mem2,0,wkMemSize);
		memset(mem3,0,wkMemSize);
		memset(pw,0,65536 * sizeof (double));

		pIn = (SSIZE *)mem1;
		for (n = 0;n < 12;n++) {
			// FFT 初期設定
			for (i = 0;i < fftSize;i++) {
				fftw_in[i][0] = pIn[((fftSize / 2) * n) + i];
				fftw_in[i][1] = 0;
			}
			// 窓関数
			for (i = 0;i < (fftSize - 1) / 2;i++) {
				fftw_in[i][0] = fftw_in[i][0] * (2.0 * i / (double)fftSize);
			}
			for (i = (fftSize - 1) / 2;i < fftSize;i++) {
				fftw_in[i][0] = fftw_in[i][0] * (2.0 - 2.0 * i / (double)fftSize);
			}

			// FFT
			fftw_execute(fftw_p);

			// 元信号のパワーを累積する
			for (i = 1;i < fftSize / 2;i++) {
				p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
				if (p != 0) {
					p = sqrt(p);
				}
				pw[i] += p;
			}
		}
#if 1
		for (i = 0;i < fftSize / 2;i++) {
			pw[i] /= 12;
		}
#endif
		for (i = 0,h = 0;h < hfc - 100;i++,h += 100) {
			// 100hz 範囲のパワーを調べる
			lowIndex  = ((double)fftSize / outSampleR) * h;
			highIndex = ((double)fftSize / outSampleR) * (h + 100);
			refPw[i] = 0;
			for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
				refPw[i] += pw[j];
			}
			if (nn > 0) {
				refPw[i] /= nn;
			}
		}
		for (n = 0;n < 3;n++) {
			// FFT 初期設定
			for (i = 0;i < fftSize;i++) {
				fftw_in[i][0] = pIn[((fftSize / 2) * n) + i];
				fftw_in[i][1] = 0;
			}
			// 窓関数
			for (i = 0;i < (fftSize - 1) / 2;i++) {
				fftw_in[i][0] = fftw_in[i][0] * (2.0 * i / (double)fftSize);
			}
			for (i = (fftSize - 1) / 2;i < fftSize;i++) {
				fftw_in[i][0] = fftw_in[i][0] * (2.0 - 2.0 * i / (double)fftSize);
			}

			// FFT
			fftw_execute(fftw_p);
			// 閾値より大きい音はカットする
			ignore_flag = 0;
			for (i = 1;i < fftSize / 2;i++) {
				hz = ((double)outSampleR / fftSize) * i;
				if (pw[i] > (double)0x100000) {
					if (ignore_flag == 0 && hz >= 100 && hz >= param->nr) {
						idx = hz / 100;
						if (pw[i] >= refPw[idx] * cutLV[param->nrLV]) {
							fftw_out[i][0] = 0;
							fftw_out[i][1] = 0;
						} else {
							fftw_out[i][0] *= 0.45;
							fftw_out[i][1] *= 0.45;
						}
					} else {
						fftw_out[i][0] = 0;
						fftw_out[i][1] = 0;
					}
				} else {
					fftw_out[i][0] = 0;
					fftw_out[i][1] = 0;
				}
			}

			// 半分のデータを復元
			for (i = 1;i < fftSize / 2;i++) {
				fftw_out[fftSize - i][0] = fftw_out[i][0];
				fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
			}
			fftw_out[0][0] = 0;
			fftw_out[0][1] = 0;

			// invert FFT
			fftw_execute(fftw_ip);

			// 出力
			pOut = (SSIZE *)mem2;
			for (i = 0;i < fftSize;i++) {
				pOut[((fftSize / 2) * n) + i] += fftw_in[i][0] / fftSize;
			}
		}
		if (startInSample + fftSize / 2 >= 0) {
			outTempFile(fp_w,mem2 + fftSize / 2,fftSize,param);
		}
	}
	free(mem1);
	free(mem2);
	free(mem3);
	free(pw);
	fftw_destroy_plan(fftw_p);
	fftw_destroy_plan(fftw_ip);
	fftw_free(fftw_in);
	fftw_free(fftw_out);

}
//---------------------------------------------------------------------------
// Function   : spAnalyze
// Description: スピーカーの周波数特性に応じて調整値パラメーターを出力する
// ---
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	param		:変換パラメータ
//
void spAnalyze(SSIZE inSample,FIO *fp_r,PARAM_INFO *param)
{
	SSIZE *mem1,*mem0;
	long outSampleR;
	long inSampleR;
	long wkMemSize;
	long fftSize,i,n,hz;
	
	
	long validIndex,adjIndex,adjWidth;
	double percent,per;
	double div,step;
	long cnt;
	double p;
	static double adjData[192000];
	SSIZE *pIn;
	SSIZE startInSample,nSample;
	fftw_complex *fftw_in;
	fftw_plan fftw_p;
	double *adjFrom,*adjTo,*adjNx;
	FILE *ofp;
	outSampleR = param->outSampleR;
	inSampleR = param->inSampleR;

	fio_rewind(fp_r);

	fftSize = param->outSampleR * 2;

	wkMemSize = (fftSize * 2) * sizeof (SSIZE);

	validIndex = ((double)fftSize / param->outSampleR) * (inSampleR / 2);
	adjIndex   = ((double)fftSize / param->outSampleR) * (2000);
	adjWidth   = ((double)fftSize / param->outSampleR) * (1000);
	if (validIndex < adjIndex) {
		adjIndex = validIndex;
	}

	adjFrom = (double *)malloc(sizeof (double) * fftSize);
	adjTo	= (double *)malloc(sizeof (double) * fftSize);
	adjNx	= (double *)malloc(sizeof (double) * fftSize);
	if (adjFrom == NULL || adjTo == NULL || adjNx == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	for (i = 0;i < fftSize;i++) {
		adjFrom[i] = 0;
		adjTo[i] = 0;
		adjNx[i] = 0;
	}
	
	mem1 = (SSIZE *)malloc(wkMemSize);
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_p = fftw_plan_dft_1d(fftSize,fftw_in,fftw_in,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	per = -1;
	cnt = 0;
	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + fftSize + (fftSize / 2);startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
			percent = ((double)startInSample / inSample);
			percent *= 100;
			percent = (double)((int)percent);
			if (percent != per) {
				fprintf(stdout,"%d%%\n",(int)percent);
				fflush(stdout);
			}
			per = percent;
//			Sleep(1);
		}

		memset(mem1,0,wkMemSize);

		if (startInSample >= 0 && startInSample + (fftSize * 2) < inSample) {
			nSample = fftSize * 2;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 2;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}

		pIn = (SSIZE *)mem1;
		for (n = 0;n < 3;n++) {
			// FFT 初期設定
			copyToFFTW(fftw_in,&pIn[((fftSize / 2) * n)],fftSize);

			// 窓関数
			windowFFTW(fftw_in,fftSize);

			// FFT
			fftw_execute(fftw_p);
			
			for (i = 1;i < validIndex;i++) {
				adjTo[i] = fftw_in[i][0] * fftw_in[i][0] + fftw_in[i][1] * fftw_in[i][1];
				if (adjTo[i] != 0) {
					adjTo[i] = sqrt(adjTo[i]);
				}
			}
			
			// 128 サイズの移動平均
			for (i = 1;i + 128 < validIndex;i++) {
				p = 0;
				for (n = 0;n < 128;n++) {
					p += adjTo[i + n];
				}
				if (p > 0) {
					p /= 128;
				}
				adjFrom[i] = p;
			}
			for (;i < fftSize / 2;i++) {
				adjFrom[i] = p;
			}

			p = 0;
			for (i = 1;i < 101;i++) {
				p += adjFrom[i];
			}
			p /= 100;
			
			adjTo[0] = 0;
			for (i = 1;i < fftSize / 2;i++) {
				adjTo[i] = p;
			}
			for (i = 1;i < fftSize / 2;i++) {
				if (adjFrom[i] != 0) {
					adjNx[i] += (adjTo[i] / adjFrom[i]);
				}
			}
			cnt++;
		}
	}
	if (cnt > 0) {
		for (i = 0;i < fftSize / 2;i++) {
			adjNx[i] /= cnt;
		}
	}
	
	p = 0;
	for (i = adjIndex,n = 0;n < adjWidth;adjIndex++,n++) {
		p += adjNx[i];
	}
	if (n > 0) {
		p /= n;
	}
	for (i = 1;i < adjIndex;i++) {
		adjNx[i] = p;
	}
	div = 1.0;
	step = 1.0 / 192000;
	for (i = adjIndex;i < fftSize / 2;i++) {
		adjNx[i] *= div;
		if (div - step > 0) {
			div -= step;
		} else {
			div = 0.08;
		}
	}

	if (param->r1_flag == 1) {
		unlink(param->sp_path);
		ofp = fopen(param->sp_path,"w");
		if (ofp) {
			for (i = 0;i < fftSize / 2;i++) {
				hz = ((double)param->outSampleR / fftSize) * i;
				adjData[hz] = adjNx[i];
			}
			for (i = 1;i < 192000;i++) {
				fprintf(ofp,"%lf\n",adjData[i]);
			}
			fclose(ofp);
		}
	}
	free(adjFrom);
	free(adjTo);
	free(adjNx);

	free(mem1);

	fftw_destroy_plan(fftw_p);
	fftw_free(fftw_in);

}

//---------------------------------------------------------------------------
// Function   : outTempFile
// Description: データをテンポラリファイルへ出力する
// ---
//	fp_w	:出力ファイル
//	in		:データのアドレス
//	size	:データー数
//	param	:パラメーター
//
void outTempFile(FIO *fp_w,void *in,SSIZE size,PARAM_INFO *param)
{
	fio_size r;

	r = fio_write(in,1,size * sizeof (SSIZE),fp_w);
	if (r != size * sizeof (SSIZE)) {
		param->err = STATUS_FILE_WRITE_ERR;
	}
}
//---------------------------------------------------------------------------
// Function   : normalNoise
// Description: 正規乱数生成
//
double normalNoise(void)
/*
 * 正規乱数
 */
{
	double x1,x2;

	x1 = (double)rand() / RAND_MAX;
	x1 = 0.99999 * x1 + 0.00001;
	x2 = (double)rand() / RAND_MAX;
	return sqrt(-log(x1)) * cos(2.0 * 3.1415926 * x2) / 3.5;
}
//---------------------------------------------------------------------------
// Function   : adjPinkFilter
// Description: 1/f 特性にするフィルター
// ---
//	mode	  :モード(0,1,2,3)
//	fftSizeOut:FFT数
//	fftw_out2 :FFTW OUT 変数
//	param	  :変換パラメータ
//
void adjPinkFilter(int mode,long fftSizeOut,fftw_complex *fftw_out2,PARAM_INFO *param)
{
	long i;
	long startIdx,endIdx;
	double hz,div,step;
	
	long cutOff;
	long outSampleR;

	outSampleR = param->outSampleR;
	if (mode == 4 && param->lpf != -1) {
		startIdx = ((double)fftSizeOut / outSampleR) * param->lpf;
		endIdx	 = ((double)fftSizeOut / outSampleR) * (param->lpf + 20000);
		step = 1.00 / ((endIdx - startIdx) * 1.55);

		div = 1;
		for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
			hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
			if (div - step > 0) {
		   		div -= step;
			} else {
		   		div = 0.01;
			}
		}
		for (;i < fftSizeOut / 2;i++) {
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
		}
		div = 1;
		for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
			hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
			if (div - step > 0) {
		   		div -= step;
			} else {
		   		div = 0.01;
			}
		}
		for (;i < fftSizeOut / 2;i++) {
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
		}
		return;
	}

	if (mode == 1) {
		// 1/f 特性にするフィルター(hfa1)
		for (i = 1;i < fftSizeOut / 2;i++) {
			hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
			if (hz > 0) {
				fftw_out2[i][0] /= hz;
				fftw_out2[i][1] /= hz;
			}
		}
	}
	if (mode != 3) {
		// hfa1、hfa2、hfa3用の高域補間時の周波数調整
		if (param->hfa != 0 && param->hfc >= 8000 && param->hfc <= 23000) {
//			if (param->hfc > 13000) {
				startIdx = ((double)fftSizeOut / outSampleR) * 13000;
				endIdx	 = ((double)fftSizeOut / outSampleR) * 27000;
				endIdx	 = ((double)fftSizeOut / outSampleR) * 24000;
				step = 1.00 / ((endIdx - startIdx) * 1.50);
				step = 1.00 / ((endIdx - startIdx) * 1.46);
				step = 1.00 / ((endIdx - startIdx) * 1.57);
				step = 1.00 / ((endIdx - startIdx) * 1.61);
				step = 1.00 / ((endIdx - startIdx) * 1.71);			// 0.8.0(TEST)
				step = 1.00 / ((endIdx - startIdx) * 1.36);
				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				if (param->validOutSampleR <= 96000) {
					startIdx = ((double)fftSizeOut / outSampleR) * 16000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 36000;
					step = 1.00 / ((endIdx - startIdx) * 1.20);
					startIdx = ((double)fftSizeOut / outSampleR) * 23000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 34000;
					step = 1.00 / ((endIdx - startIdx) * 1.46);
					step = 1.00 / ((endIdx - startIdx) * 1.61);
				} else {
					startIdx = ((double)fftSizeOut / outSampleR) * 23000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 34000;
					if (param->overSamp < 1) {
						startIdx = ((double)fftSizeOut / outSampleR) * 24000;	// 0.8.0(TEST)
						endIdx	 = ((double)fftSizeOut / outSampleR) * 64000;	// 0.8.0(TEST)
						step = 1.00 / ((endIdx - startIdx) * 1.66);
						step = 1.00 / ((endIdx - startIdx) * 1.50);
						step = 1.00 / ((endIdx - startIdx) * 1.53);				// 0.8.0(TEST)
					} else {
						startIdx = ((double)fftSizeOut / outSampleR) * 24000;	// 0.8.0(TEST2)
						endIdx	 = ((double)fftSizeOut / outSampleR) * 48000;
						step = 1.00 / ((endIdx - startIdx) * 1.48);				// 0.8.0(TEST2)
					}
				}
				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				if (param->validOutSampleR <= 96000) {
					startIdx = ((double)fftSizeOut / outSampleR) * 22000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 68000;
//					step = 1.00 / ((endIdx - startIdx) * 1.2);
					step = 1.00 / ((endIdx - startIdx) * 1.25);
					startIdx = ((double)fftSizeOut / outSampleR) * 31000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 47000;
					step = 1.00 / ((endIdx - startIdx) * 1.61);
					step = 1.00 / ((endIdx - startIdx) * 1.68);
					step = 1.00 / ((endIdx - startIdx) * 1.71);
				} else {
					startIdx = ((double)fftSizeOut / outSampleR) * 29000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 57000;
					startIdx = ((double)fftSizeOut / outSampleR) * 31000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 54000;
					if (param->overSamp < 1) {
						startIdx = ((double)fftSizeOut / outSampleR) * 66000;		// 0.8.0(TEST)
						endIdx	 = ((double)fftSizeOut / outSampleR) * 120000;		// 0.8.0(TEST)
						step = 1.00 / ((endIdx - startIdx) * 1.61);
						step = 1.00 / ((endIdx - startIdx) * 1.21);
						step = 1.00 / ((endIdx - startIdx) * 1.53);					// 0.8.0(TEST)
					} else {
						startIdx = ((double)fftSizeOut / outSampleR) * 48000;		// 0.8.0(TEST2)
						endIdx	 = ((double)fftSizeOut / outSampleR) * 120000;		// 0.8.0(TEST2)
						step = 1.00 / ((endIdx - startIdx) * 1.65);					// 0.8.0(TEST2)
					}
				}
				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
#if 0
				if (param->validOutSampleR <= 96000) {
					startIdx = ((double)fftSizeOut / outSampleR) * 43000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 68000;
					step = 1.00 / ((endIdx - startIdx) * 1.71);
				} else {
					startIdx = ((double)fftSizeOut / outSampleR) * 43000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 78000;
					step = 1.00 / ((endIdx - startIdx) * 1.81);
				}
				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
#endif
				if (param->validOutSampleR <= 96000) {
					startIdx = ((double)fftSizeOut / outSampleR) * 30000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 120000;
//					step = 1.00 / ((endIdx - startIdx) * 1.2);
					step = 1.00 / ((endIdx - startIdx) * 1.25);
					startIdx = ((double)fftSizeOut / outSampleR) * 55000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 88000;
					step = 1.00 / ((endIdx - startIdx) * 1.81);
				} else {
					startIdx = ((double)fftSizeOut / outSampleR) * 70000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 120000;
//					step = 1.00 / ((endIdx - startIdx) * 1.2);
					step = 1.00 / ((endIdx - startIdx) * 1.25);
					startIdx = ((double)fftSizeOut / outSampleR) * 71000;
					endIdx	 = ((double)fftSizeOut / outSampleR) * 120000;
					if (param->overSamp < 1) {
						startIdx = ((double)fftSizeOut / outSampleR) * 77000;			// 0.8.0(TEST)
						endIdx	 = ((double)fftSizeOut / outSampleR) * 120000;			// 0.8.0(TEST)
						step = 1.00 / ((endIdx - startIdx) * 1.21);
						step = 1.00 / ((endIdx - startIdx) * 1.71);						// 0.8.0(TEST)
					} else {
						startIdx = ((double)fftSizeOut / outSampleR) * 77000;			// 0.8.0(TEST2)
						endIdx	 = ((double)fftSizeOut / outSampleR) * 95000;			// 0.8.0(TEST2)
						step = 1.00 / ((endIdx - startIdx) * 1.21);						// 0.8.0(TEST2)
					}
				}
				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
#if 0
			} else {
				startIdx = ((double)fftSizeOut / outSampleR) * 10000;
				endIdx	 = ((double)fftSizeOut / outSampleR) * 20000;

				step = 1.00 / ((endIdx - startIdx) * 1.68);

				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}

				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
			}
#endif
		}

		if (mode == 2) {
			// 独自のローパスフィルター
			cutOff = 40000;
			if (param->lpf > 1 && cutOff > param->lpf) {
				cutOff = param->lpf;
			}
			if ((outSampleR / 2) >= cutOff) {
				startIdx  = ((double)fftSizeOut / outSampleR) * cutOff;
				if (param->overSamp == 0) {
					endIdx = fftSizeOut / 2;
				} else {
					endIdx = fftSizeOut / 4;
				}
				step = 1.00 / ((endIdx - startIdx) * 1.30);

				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
					if (div - step > 0) {
						div -= step;
					} else {
						div = 0.01;
					}
				}
			}
		}
	} else {
		// デエンファシス用の処理
		if (param->deEmphasis == 1) {
			startIdx = ((double)fftSizeOut / outSampleR) * 3180;
			endIdx	 = ((double)fftSizeOut / outSampleR) * 10600;
		} else {
			startIdx = ((double)fftSizeOut / outSampleR) * 2100;
			endIdx	 = ((double)fftSizeOut / outSampleR) * 9520;
		}
		step = 1.00 / ((endIdx - startIdx) * 1.75);

		div = 1;
		for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
			hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
			if (div - step > 0) {
		   		div -= step;
			} else {
		   		div = 0.01;
			}
		}
		for (;i < fftSizeOut / 2;i++) {
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
		}
	}

	if (param->overSamp == 0) {
		for (i = (fftSizeOut / 2) - 5;i < fftSizeOut / 2;i++) {
			fftw_out2[i][0] = 0;
			fftw_out2[i][1] = 0;
		}
	} else {
		for (i = (fftSizeOut / 4) - 5;i < fftSizeOut / 4;i++) {
			fftw_out2[i][0] = 0;
			fftw_out2[i][1] = 0;
		}
	}
}
//---------------------------------------------------------------------------
// Function   : merageTempFile
// Description: 出力結果のファイルをマージする
// ---
//	type	 :マージの種類
//	normFlag :ノーマライズ用変数更新フラグ
//	fp_r 	 :入力ファイル1
//	fp_r2	 :入力ファイル2
//	fp_w	 :出力ファイル
//	inSample :サンプル数
//	param	 :パラメーター
//
void merageTempFile(char type,int normFlag,FIO *fp_r,FIO *fp_r2,FIO *fp_w,SSIZE inSample,PARAM_INFO *param)
{
	
	SSIZE min,max;
	SSIZE maxLv,maxLv2;
	SSIZE remainSample;
	SSIZE ns;
	long i;
	fio_size remain1,remain2;
	fio_size wr_n;

	fio_rewind(fp_r);

	if (fp_r2 != NULL) {
		fio_rewind(fp_r2);
	}

	if (fp_w != NULL) {
		fio_rewind(fp_w);
	}

	min = max = 0;
	ns	= 0;
	maxLv = 0;
	maxLv2 = 0;
	remainSample = inSample;

	do {
		if (type == '+') {
//			Sleep(1);
			remain1 = fio_read(diskBuffer,sizeof (SSIZE),1 * 1024 * 1024,fp_r);
			if (fp_r2 != NULL) {
				remain2 = fio_read(diskBuffer2,sizeof (SSIZE),1 * 1024 * 1024,fp_r2);
			}
			if (remain1 == 0 || remain2 == 0) {
				break;
			}
			for (i = 0;i < remain1;i++) {
				if (diskBuffer[i] != 0) {
					diskBuffer[i] += diskBuffer2[i];
					if (diskBuffer[i] < min) {
						min = diskBuffer[i];
					}
					if (diskBuffer[i] > max) {
						max = diskBuffer[i];
					}
					if (diskBuffer[i] >> 40 > 0) {
						maxLv2 += diskBuffer[i] >> 40;
						ns++;
						if (maxLv2 >= 0x1000000000000) {
							maxLv2 /= ns;
							if (maxLv > 0) {
								maxLv = (maxLv + maxLv2) / 2;
							} else {
								maxLv = maxLv2;
							}
							maxLv2 = 0;
							ns = 0;
						}
					}
				}
			}
			remain1 =  remain1 < remainSample ? remain1 : remainSample;
			if (fp_w != NULL) {
//				fio_seek(fp_w,pos,SEEK_SET);
				wr_n = fio_write(diskBuffer,sizeof (SSIZE),remain1,fp_w);
				if (wr_n != remain1) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
			}
			remainSample -= remain1;
		} else if (type == '-') {
//			Sleep(1);
			remain1 = fio_read(diskBuffer,sizeof (SSIZE),1 * 1024 * 1024,fp_r);
			if (fp_r2 != NULL) {
				remain2 = fio_read(diskBuffer2,sizeof (SSIZE),1 * 1024 * 1024,fp_r2);
			}
			if (remain1 == 0) {
				break;
			}
			for (i = 0;i < remain1;i++) {
				if (diskBuffer[i] != 0) {
					diskBuffer[i] -= diskBuffer2[i];
					if (diskBuffer[i] < min) {
						min = diskBuffer[i];
					}
					if (diskBuffer[i] > max) {
						max = diskBuffer[i];
					}
					if (diskBuffer[i] >> 40 > 0) {
						maxLv2 += diskBuffer[i] >> 40;
						ns++;
						if (maxLv2 >= 0x1000000000000) {
							maxLv2 /= ns;
							if (maxLv > 0) {
								maxLv = (maxLv + maxLv2) / 2;
							} else {
								maxLv = maxLv2;
							}
							maxLv2 = 0;
							ns = 0;
						}
					}
				}
			}
			remain1 =  remain1 < remainSample ? remain1 : remainSample;
			if (fp_w != NULL) {
//				fio_seek(fp_w,pos,SEEK_SET);
				wr_n = fio_write(diskBuffer,sizeof (SSIZE),remain1,fp_w);
				if (wr_n != remain1) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
			}
			remainSample -= remain1;
		} else if (type == ' ') {
			//Sleep(1);
			remain1 = fio_read(diskBuffer,sizeof (SSIZE),1 * 1024 * 1024,fp_r);
			if (remain1 == 0) {
				break;
			}
			for (i = 0;i < remain1;i++) {
				if (diskBuffer[i] < min) {
					min = diskBuffer[i];
				}
				if (diskBuffer[i] > max) {
					max = diskBuffer[i];
				}
				if (diskBuffer[i] >> 40 > 0) {
					maxLv2 += diskBuffer[i] >> 40;
					ns++;
					if (maxLv2 >= 0x1000000000000) {
						maxLv2 /= ns;
						if (maxLv > 0) {
							maxLv = (maxLv + maxLv2) / 2;
						} else {
							maxLv = maxLv2;
						}
						ns = 0;
					}
				}
			}
			remain1 =  remain1 < remainSample ? remain1 : remainSample;
			if (fp_w != NULL) {
//				fio_seek(fp_w,pos,SEEK_SET);
				wr_n = fio_write(diskBuffer,sizeof (SSIZE),remain1,fp_w);
				if (wr_n != remain1) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
			}
			remainSample -= remain1;
		}
	} while (remain1 == 1 * 1024 * 1024L && remainSample > 0);

	if (remainSample > 0 && param->err == STATUS_SUCCESS) {
		param->err = STATUS_FILE_READ_ERR;param->errLine = __LINE__;
	}

	if (param->err) {
		param->err = STATUS_FILE_READ_ERR;
		// エラー終了
		return;
	}
	
	if (fp_w != NULL) {
		fio_flush(fp_w);
	}
	if (normFlag == 1) {
		if (max > NormInfo.max) {
			NormInfo.max = max;
		}
		if (min < NormInfo.min) {
			NormInfo.min = min;
		}
		if (ns > 0) {
			maxLv2 /= ns;
		}
		if (maxLv > 0) {
			maxLv = (maxLv + maxLv2) / 2;
		} else {
			maxLv = maxLv2;
		}
		NormInfo.avg = maxLv;
	}
}
#if 0
//---------------------------------------------------------------------------
// Function   : bpFilter
// Description: 指定周波数をカットする
// ---
//	lfc		 	:低域のカットオフ周波数(この周波数以下の領域をカットする)
//	hfc		 	:高域のカットオフ周波数(この周波数以上の領域をカットする)
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	fp_w		:出力ファイル
//	param		:変換パラメータ
//
int bpFilter(long lfc,long hfc,DWORD inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3,*mem4;
	long outSampleR;
	long wkMemSize;
	long fftSize,i;
	SSIZE *pIn[3],*pOut[3];
	SSIZE startInSample,nSample;
	fftw_complex *fftw_in[3],*fftw_out[3];
	fftw_plan fftw_p[3],fftw_ip[3];

	outSampleR = param->outSampleR;
	fio_rewind(fp_r);
	fio_rewind(fp_w);

	fftSize = param->outSampleR * 2;

	wkMemSize = (fftSize * 2) * sizeof (SSIZE);

	mem1 = (SSIZE *)al_malloc(wkMemSize);
	if (mem1 == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	mem2 = (SSIZE *)al_malloc(wkMemSize);
	if (mem2 == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	mem3 = (SSIZE *)al_malloc(wkMemSize);
	if (mem3 == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	mem4 = (SSIZE *)al_malloc(wkMemSize);
	if (mem4 == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}

	fftw_in[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[0] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_in[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[1] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_in[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[2] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}

	fftw_out[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[0] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_out[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[1] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_out[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[2] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}

	fftw_p[0] = fftw_plan_dft_1d(fftSize,fftw_in[0],fftw_out[0],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[0] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_p[1] = fftw_plan_dft_1d(fftSize,fftw_in[1],fftw_out[1],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[1] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_p[2] = fftw_plan_dft_1d(fftSize,fftw_in[2],fftw_out[2],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[2] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}

	fftw_ip[0] = fftw_plan_dft_1d(fftSize,fftw_out[0],fftw_in[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[0] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_ip[1] = fftw_plan_dft_1d(fftSize,fftw_out[1],fftw_in[1],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[1] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_ip[2] = fftw_plan_dft_1d(fftSize,fftw_out[2],fftw_in[2],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[2] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}


	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + fftSize + (fftSize / 2);startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
//			Sleep(0);
		}

		memset(mem1,0,wkMemSize);

		if (startInSample >= 0 && startInSample + (fftSize * 2) < inSample) {
			nSample = fftSize * 2;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 2;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}

		pIn[0]	= &mem1[((fftSize / 2) * 0)];
		pOut[0] = &mem2[((fftSize / 2) * 0)];
		pIn[1]	= &mem1[((fftSize / 2) * 1)];
		pOut[1] = &mem3[((fftSize / 2) * 1)];
		pIn[2]	= &mem1[((fftSize / 2) * 2)];
		pOut[2] = &mem4[((fftSize / 2) * 2)];

		memset(mem2,0,wkMemSize);
		memset(mem3,0,wkMemSize);
		memset(mem4,0,wkMemSize);

		#pragma omp parallel
		{
			#pragma omp sections
			{
				#pragma omp section
				{
					// 1
					bpFilterSub(pIn[0],pOut[0],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],lfc,hfc,param);
				}
				#pragma omp section
				{
					// 2
					bpFilterSub(pIn[1],pOut[1],fftw_in[1],fftw_out[1],fftw_p[1],fftw_ip[1],lfc,hfc,param);
				}
				#pragma omp section
				{
					// 3
					bpFilterSub(pIn[2],pOut[2],fftw_in[2],fftw_out[2],fftw_p[2],fftw_ip[2],lfc,hfc,param);
				}
			}
			#pragma omp for
			for (i = 0;i < nSample;i++) {
				mem2[i] += mem3[i] + mem4[i];
			}
		}

		if (startInSample + fftSize / 2 >= 0) {
			fio_seek(fp_w,(startInSample + (fftSize / 2)) * sizeof (SSIZE),SEEK_SET);
			outTempFile(fp_w,mem2 + fftSize / 2,fftSize,param);
			if (param->err) {
				break;
			}
		}
	}

	fio_flush(fp_w);

	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	al_free(mem4);

	fftw_destroy_plan(fftw_p[0]);
	fftw_destroy_plan(fftw_p[1]);
	fftw_destroy_plan(fftw_p[2]);

	fftw_destroy_plan(fftw_ip[0]);
	fftw_destroy_plan(fftw_ip[1]);
	fftw_destroy_plan(fftw_ip[2]);

	fftw_free(fftw_in[0]);
	fftw_free(fftw_in[1]);
	fftw_free(fftw_in[2]);

	fftw_free(fftw_out[0]);
	fftw_free(fftw_out[1]);
	fftw_free(fftw_out[2]);

	return STATUS_SUCCESS;
}
//---------------------------------------------------------------------------
// Function   : bpFilterSub
// Description: 指定周波数をカットする
// ---
//	param		:変換パラメータ
//
void bpFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,long lfc,long hfc,PARAM_INFO *param)
{
	long fftSize;
	long lowIndex,highIndex;
	long outSampleR;
	long i;

	fftSize = param->outSampleR * 2;
	outSampleR = param->outSampleR;
	
	// FFT 初期設定
	copyToFFTW(fftw_in,pIn,fftSize);

	// 窓関数
	windowFFTW(fftw_in,fftSize);

	// FFT
	fftw_execute(fftw_p);

	// 元信号の高域のパワーを調べる
	if (lfc != -1) {
		lowIndex = ((double)fftSize / outSampleR) * (lfc);
		
		// 低域カット
		for (i = 1;i < lowIndex;i++) {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
	}
	if (hfc != -1) {
		highIndex = ((double)fftSize / outSampleR) * (hfc);
		for (i = highIndex;i < fftSize / 2;i++) {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
	}
	// 半分のデータを復元
	for (i = 1;i < fftSize / 2;i++) {
		fftw_out[fftSize - i][0] = fftw_out[i][0];
		fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
	}

	// invert FFT
	fftw_execute(fftw_ip);

	// 出力
	for (i = 0;i < fftSize;i++) {
		pOut[i] += fftw_in[i][0] / fftSize;
	}
}
#endif
static int chkAbort(PARAM_INFO *param,int percent,int percent_diff)
{
	FILE *fp;
	int diff = param->abort_percent - percent;
	
	if (diff == 0 || diff > percent_diff || diff < (percent_diff * -1) || percent == 0 || percent == 100) {
		param->abort_percent = percent;
		if (param->abort_filename[0] != '\0') {
			fp = fopen(param->abort_filename,"r");
			if (fp) {
				fclose(fp);
				return 1;
			}
		}
	}
	return 0;
}

