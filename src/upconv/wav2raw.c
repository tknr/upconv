/****************************************************************************/
/* wav2raw (C) 2011-2019 By 59414d41										*/
/* wavファイルを64bitのrawファイルへ変換する								*/
/* rawファイルをwavファイルへ変換する										*/
/* dsfファイルを64bitのrawファイルへ変換する								*/
/* rawファイルをdsfファイルにする											*/
/*																			*/
/****************************************************************************/

/*--- Log ------------------------------------------------------------------
 * Ver 0.10 <09/07/15> - upconvから分離
 * Ver 0.21 <09/10/26> - ちょっと修正
 * Ver 0.30 <09/11/01> - パラメータファイルの採用
 * Ver 0.31 <09/11/16> - エラー情報をファイルへ出力するようにした
 * Ver 0.50 <10/11/02> - 処理修正
 * Ver 0.70 <11/07/24> - コンパイラをmingwに変更
 *						 大きなファイルに対応
 *						 split処理をやめたことによる修正
 * Ver 0.80 <12/02/11> - fileioを使用するように修正
 *						 マルチチャンネルに対応
 * ver 0.99 <18/10/25> - raw2wav,dsf2raw とマージ
   ver 1.20 <19/10/12> - upconv.c から呼び出すように修正
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <time.h>
#include "upconv.h"
#include "fftw3.h"
#include "./fileio.h"
#include "./../PLG_AUDIO_IO/PLG_AudioIO.h"

#if 1
#define	PRINT_LOG(s)	do {																	\
							FILE *log;															\
							log = fopen("d:\\wav2raw.log","a");									\
							if (log) {															\
								fprintf(log,"%s [%d] %s\n",__FUNCTION__,__LINE__,s);			\
								fclose(log);													\
							}																	\
						} while (0)
#else
#define	PRINT_LOG(s)	//
#endif

// サンプルを処理するデータ型
#define SSIZE	signed long long int
#define UI64	unsigned long long int

typedef struct {
	SSIZE	minR1;
	SSIZE	maxR1;
	SSIZE	minR2;
	SSIZE	maxR2;
	SSIZE	minR3;
	SSIZE	maxR3;
	SSIZE	minR4;
	SSIZE	maxR4;
	SSIZE	minR5;
	SSIZE	maxR5;
	SSIZE	minR6;
	SSIZE	maxR6;
	SSIZE	maxLv1;
	SSIZE	maxLv2;
	SSIZE	maxLv3;
	SSIZE	maxLv4;
	SSIZE	maxLv5;
	SSIZE	maxLv6;
	SSIZE	tmpLv1;
	SSIZE	tmpLv2;
	SSIZE	tmpLv3;
	SSIZE	tmpLv4;
	SSIZE	tmpLv5;
	SSIZE	tmpLv6;
	SSIZE	cntLv1;
	SSIZE	cntLv2;
	SSIZE	cntLv3;
	SSIZE	cntLv4;
	SSIZE	cntLv5;
	SSIZE	cntLv6;
} NORM_INFO;

typedef struct {
	long sampling;
	long bitwidth;
	long mode;
	long norm;
	long norm_option;
	long encorder;
	long split;
	long downmix;
	long	dsd_fmt;
	int ch;
	int chC;
	int chS;
	int chLFE;
	int bwf;
	int raw;
	int rf64;
	int abe;
	int hfa;
	int flac;
	int thread;
	int fio;
	int	errLine;
	double	nx1;
	double	nx2;
	double	nx3;
	double	nx4;
	double	nx5;
	double	nx6;
	SSIZE	sum1;
	SSIZE	sum2;
	SSIZE	sum3;
	SSIZE	sum4;
	SSIZE	sum5;
	SSIZE	sum6;
	SSIZE	sum1_2nd;
	SSIZE	sum2_2nd;
	SSIZE	sum3_2nd;
	SSIZE	sum4_2nd;
	SSIZE	sum5_2nd;
	SSIZE	sum6_2nd;
	SSIZE	sum1_3rd;
	SSIZE	sum2_3rd;
	SSIZE	sum3_3rd;
	SSIZE	sum4_3rd;
	SSIZE	sum5_3rd;
	SSIZE	sum6_3rd;
	SSIZE	old_s1;
	SSIZE	old_s2;
	SSIZE	old_s3;
	SSIZE	old_s4;
	SSIZE	old_s5;
	SSIZE	old_s6;
	SSIZE	s_b1,s_a1;
	SSIZE	s_b2,s_a2;
	SSIZE	s_b3,s_a3;
	SSIZE	s_b4,s_a4;
	SSIZE	s_b5,s_a5;
	SSIZE	s_b6,s_a6;
	int ditherLv;
	char fromfile[_MAX_PATH];
	char tofile[_MAX_PATH];
	char opt_flac[512];
	char opt_wavpack[512];
	char opt_mp3[512];
	FIO *fp_w[6];
	FIO *fp_temp[2];
} PARAM_INFO;

typedef struct {
	UI64	sampling;
	UI64	data_offset;
	UI64	n_sample;
	int		channel;
	int		err;
	int		errLine;
	int		decode;
	long	thread;
	long	inSampleR;
	long	outSampleR;
	long	enable_hfc;
	long	lfc;
	long	hfc;
	long	src_flag;
	long	dsf_mode;
	long	fio;
	long	dsd_fmt;
	SSIZE	l_min,l_max;
	SSIZE	r_min,r_max;
	SSIZE	l_maxLv,r_maxLv;
	SSIZE	l_tmpLv,r_tmpLv;
	SSIZE	l_cntLv,r_cntLv;
	char	*workpath;
	char	*argv4;
} PARAM_INFO2;

// BWF の link チャンク対応
char link_start[]="<LINK>\r\n";
char link_file[]=						\
	"\t<FILE type=\"%s\">\r\n"				\
	"\t\t<FILENUMBER>%d</FILENUMBER>\r\n"	\
	"\t\t<FILENAME>%s</FILENAME>\r\n"		\
	"\t</FILE>\r\n";
char link_end[]="</LINK>";

//
// DSF ファイルフォーマット仕様書を参照
#pragma pack(push, 1)
typedef struct {
	char	id[4];
	UI64	chunk_size;
	UI64	file_size;
	UI64	ptr_meta;
} DSF;

typedef struct {
	char	id[4];
	UI64	chunk_size;
	DWORD	fmt_version;
	DWORD	fmt_id;
	DWORD	channel_type;
	DWORD	channel_count;
	DWORD	sampling;
	DWORD	sample_bit_count;
	UI64	sample_count;
	DWORD	block_size;
	DWORD	reserved;
} DSF_FMT;

typedef struct {
	char	id[4];
	UI64	chunk_size;
} DSF_DATA;
#pragma pack(pop)

NORM_INFO NormInfo;
PARAM_INFO paramInfo;

int errLine;

/*--- Function Prototype ---------------------------------------------------*/
int to_raw_main(int argc, char *argv[]);
int dsf_main(int argc, char *argv[]);
int to_wav_main(int argc, char *argv[]);
extern double normalNoise(void);
extern int start_exec(int argc,char *argv[],int cpu_pri,HANDLE *ret_hStdOutRead,HANDLE *ret_process_id,HANDLE *ret_thread_id);

/*--- Wav Function Prototype ---------------------------------------------------*/
int Normalize(int *nCount,SOUNDFMT *inFmt,SOUNDFMT *outFmt,FIO *fp_r1,FIO *fp_r2,FIO *fp_r3,FIO *fp_r4,FIO *fp_r5,FIO *fp_r6,DWORD *start,DWORD nSample,PARAM_INFO *param);
int Normalize_Mx(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int Normalize_M0(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int Normalize_M1(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int Normalize_M2(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int Normalize_M3(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int UpdateBext(BROADCAST_EXT *bext,SOUNDFMT *inFmt,SOUNDFMT *outFmt,PARAM_INFO *param,long bwf_size);
int nBitTo64S(int nCh,int ch,int bit,void *in,FIO *fio,DWORD nSample);

/*--- DSF Function Prototype ---------------------------------------------------*/
int Normalize_DSD(SOUNDFMT *inFmt,SOUNDFMT *outFmt,FIO *fp_r1,FIO *fp_r2,DWORD nSample,PARAM_INFO *param);
void dsf_encode(char *in_file,char *out_file,PARAM_INFO2 *param);
void dsf_decode(char *in_file,char *out_file,PARAM_INFO2 *param);
static void fftFilter(int lr,SSIZE inSample,SSIZE outSample,FIO *fp_r,FIO *fp_w,PARAM_INFO2 *param);
static void fftFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,PARAM_INFO2 *param,int id);
void onebit2nbit(SSIZE offset,SSIZE n,SSIZE *buffer,FIO *fp_r,PARAM_INFO2 *param);
void deinterleave(UI64 inByte,FIO *fp_r,FIO *fp_w1,FIO *fp_w2,PARAM_INFO2 *param);
void nbit2onebit(SSIZE *i_buf,BYTE *o_buf,int size);
void ana_abe(SSIZE start,UI64 nSample,SSIZE *i_buffer,SSIZE *o_buffer,PARAM_INFO2 *param);
void copyToFFTW(fftw_complex *fftw,SSIZE *buf,long size);
void windowFFTW(fftw_complex *fftw,long size);
void cutFFTW(fftw_complex *fftw,long index,long size);
void *al_malloc(long size);
void *al_free(void *ptr);

double normalNoise(void);

#if 1
#define CLIP_NX(v,nx)	(v) == 0 ? ((SSIZE)0) : \
							((v) > 0 ?	\
								(((SSIZE)(v) * nx) >= (SSIZE)(0x007FFFFFFFFFFFFF) ?		(SSIZE)(0x7FFFFFFFFFFFFF)		: ((SSIZE)((v) * nx))) : \
								(((SSIZE)(v) * nx) <= ((SSIZE)0x007FFFFFFFFFFFFF * -1) ? ((SSIZE)0x7FFFFFFFFFFFFF * -1) : ((SSIZE)((v) * nx))))
#else
#define CLIP_NX(v,nx)	((SSIZE)((v) * (nx)) << 6)
#endif

#define ROUND_NBIT(b)	((SSIZE)1 << (b) - 1) 

#define CLIP_ADD(v,a)	(SSIZE)(v) + (a)

#define CLIP_MAX(v)	(v) == 0 ? ((SSIZE)0) : \
							((v) > 0 ?	\
								((SSIZE)(v) >= (SSIZE)(0x007FFFFFFFFFFFFF) ?		(SSIZE)(0x7FFFFFFFFFFFFFFF)		: ((SSIZE)((v) << 8))) : \
								(((SSIZE)(v)) <= ((SSIZE)0x007FFFFFFFFFFFFF * -1) ? ((SSIZE)0x7FFFFFFFFFFFFFFF * -1) : ((SSIZE)((v) << 8))))

#define CLIP_MAX_N(v,b)	(v) == 0 ? ((SSIZE)0) : \
							((v) > 0 ?	\
								((SSIZE)(v) >= ((((SSIZE)1 << (b)) - 1) >> 1) ?		((((SSIZE)1 << (b)) - 1) >> 1)		: (SSIZE)(v)) : \
								(((SSIZE)(v)) <= (((((SSIZE)1 << (b)) - 1) * -1) >> 1) ? (((((SSIZE)1 << (b)) - 1) * -1) >> 1) : ((SSIZE)(v))))

//---------------------------------------------------------------------------
// Function   : main
// Description: 引数を処理し変換関数を呼び出す
//
//
#if 0
int main(int argc, char *argv[])
{
	int retCode;
	int e;

	e = 1;
	if (argc == 4) {
		if (strcmp(argv[1],"-d") == 0) {
			retCode = d_main(argc,argv);
			e = 0;
		} else if (strcmp(argv[1],"-dsf") == 0) {
			retCode = dsf_main(argc,argv);
			e = 0;
		} else if (strcmp(argv[1],"-e") == 0) {
			retCode = e_main(argc,argv);
			e = 0;
		}
	}
	if (e) {
		printf(STR_COPYRIGHT);
		printf(STR_USAGE);
		exit(0);
	}
}
#endif

//---------------------------------------------------------------------------
// Function   : to_raw_main
// Description: upconv 処理をするため、raw ファイルへ変換する
//
// argv[1] Input  WAV File
// argv[2] Output WAV File
// argv[3] default parameter
// argv[4] parameter
//

int to_raw_main(int argc, char *argv[])
{
	static char workpath[_MAX_PATH];
	static char tmppath[_MAX_PATH];
	static char drive[_MAX_DRIVE];
	static char dir[_MAX_DIR];
	static char workdrive[_MAX_DRIVE];
	static char workdir[_MAX_DIR];
	static char workfname[_MAX_FNAME];
	static char workext[_MAX_EXT];
	char fname[_MAX_FNAME];
	char fname2[_MAX_FNAME];
	char ext[_MAX_EXT];
	char work[2048];
	char pparam[2048];
	char *p1,*p2;
	long temp;
	FIO  fp_r;
	FIO  fp_w1,fp_w2,fp_w3,fp_w4,fp_w5,fp_w6;
	FILE *ofp;
	FILE *fp_param;
	FILE *fp_files;
	SOUNDFMT inFmt;
	unsigned char *inBuffer;
	DWORD i,inSample,nSample;
	fio_size rs,rd;
	fio_size wr;
	fio_size max_size;
	int retCode;
	FILEINFO fileInfo;
	STARTEXEC_INFO startexec_info;
	fio_size seekPtr;
	fio_size rd_byte;
	SSIZE max,avg;
	double persent,per;
	int err1,err2,err3,err4,err5,err6;
	int thread = 1;
	int fio;
	PARAM_INFO2 param_info2;

	fio = 5;
	errLine = 0;

	pparam[0] = '\0';
	if (argc == 5) {
		do {
			memset(&NormInfo,0,sizeof (NORM_INFO));
			memset(&fileInfo,0,sizeof (FILEINFO));
			memset(&startexec_info,0,sizeof (STARTEXEC_INFO));
			memset(&param_info2,0,sizeof (PARAM_INFO2));
			param_info2.l_min = 0;
			param_info2.l_max = 0;
			param_info2.r_min = 0;
			param_info2.r_max = 0;
			param_info2.thread = 1;
			param_info2.n_sample = 0;
			param_info2.dsf_mode = 0;
			param_info2.fio = -1;
			param_info2.hfc = -1;
			param_info2.enable_hfc = 0;

			// default parameter ファイル
			fp_param = fopen(argv[3],"r");
			if (fp_param == NULL) {
				retCode = STATUS_PARAMETER_ERR;errLine = __LINE__;
				break;
			}
			
			// パラメータの読みこみ
			if (fgets(work,2047,fp_param) == NULL) {
				retCode = STATUS_PARAMETER_ERR;errLine = __LINE__;
				break;
			}
			p1 = strrchr(work,'\n');if (p1 != NULL) *p1 = '\0';
			strcat(pparam,work);
			strcat(pparam," ");
			if (strlen(argv[4]) >= 1) strcat(pparam,argv[4]);

			// tmpファイル用の作業ディレクトリ
			if (fgets(workpath,_MAX_PATH,fp_param) == NULL) {
				retCode = STATUS_PARAMETER_ERR;errLine = __LINE__;
				break;
			}
			p1 = strrchr(workpath,'\n');if (p1 != NULL) *p1 = '\0';
			if (strlen(workpath) >= 2 && workpath[strlen(workpath) - 1] != '\\') strcat(workpath,"\\");

			param_info2.workpath = strdup(workpath);
			fclose(fp_param);
			fp_param = NULL;

			_splitpath(argv[2],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,fname,"err");
			unlink(tmppath);

			_splitpath(argv[2],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,fname,"files");
			fp_files = fopen(tmppath,"w");
			if (fp_files == NULL) {
				retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
				break;
			}

			_splitpath(argv[2],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,fname,"param");
			fp_param = fopen(tmppath,"w");
			if (fp_param == NULL) {
				retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
				break;
			}
			fprintf(fp_files,"%s\n",tmppath);

			p1 = pparam;
			p2 = strchr(p1,(int)' ');

			for (;p1 != NULL;) {
				if (p2 != NULL) *p2 = '\0';

				if (sscanf(p1,"-thread:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 24) {
						thread = (int)temp;
						param_info2.thread = (int)temp;
					}
				}
				if (sscanf(p1,"-fio:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 10000) {
						fio = temp;
						param_info2.fio = temp;
						
					}
				}
				if (strcmp(p1,"-enable_hfc") == 0) {
					param_info2.enable_hfc = 1;
				}
				if (sscanf(p1,"-hfc:%ld",&temp) == 1) {
					if (temp >= 1000 && temp <= (384000 / 2)) {
						param_info2.hfc = temp;
					}
				}
				if (sscanf(p1,"-dsf:%ld",&temp) == 1) {
					if (temp >= 0 && temp <= 2) {
						param_info2.dsf_mode = (int)temp;
					}
				}
				if (p2 == NULL) {
					break;
				}
				p1 = p2 + 1;
				p2 = strchr(p1,(int)' ');
			}
			if (param_info2.enable_hfc == 0) {
				param_info2.hfc = -1;
			}

			// 音声ファイル情報を取得する
			retCode = PLG_InfoAudioData(argv[1],&inFmt,&inSample,&fileInfo);
			if (retCode != STATUS_SUCCESS) {
				retCode = STATUS_FILE_READ_ERR;errLine = __LINE__;
				break;
			}
			
			if (inFmt.channel < 1 || inFmt.channel > 6) {
				retCode = STATUS_FILE_READ_ERR;errLine = __LINE__;
				break;
			}
			
			if (strcmp(inFmt.fmt,"wav") != 0 && strcmp(inFmt.fmt,"dsf") != 0 && strcmp(inFmt.fmt,"rf64") != 0) {
				char *arg[4];
				char arg0[4096];
				char arg1[4096];
				char arg2[4096];
				char arg3[4096];
				char dec_wav[4096];
				_splitpath(argv[0],drive,dir,fname,ext);
				if (strcmp(inFmt.fmt,"flac") == 0) {
					sprintf(arg0,"%s%s%s",drive,dir,"flac.exe");
					sprintf(arg1,"|-d -o");
					sprintf(arg2,"%s.wav",argv[2]);
					sprintf(arg3,"%s",argv[1]);
				} else if (strcmp(inFmt.fmt,"wavp") == 0) {
					sprintf(arg0,"%s%s%s",drive,dir,"wvunpack.exe");
					sprintf(arg1,"-q");
					sprintf(arg2,"%s",argv[1]);
					sprintf(arg3,"%s.wav",argv[2]);
				} else if (strcmp(inFmt.fmt,"mp3") == 0) {
					sprintf(arg0,"%s%s%s",drive,dir,"lame.exe");
					sprintf(arg1,"--decode");
					sprintf(arg2,"%s",argv[1]);
					sprintf(arg3,"%s.wav",argv[2]);
				}
				strcpy(dec_wav,argv[2]);strcat(dec_wav,".wav");
				arg[0] = arg0;
				arg[1] = arg1;
				arg[2] = arg2;
				arg[3] = arg3;
				sprintf(work,"argv:[%s,%s,%s,%s]",arg0,arg1,arg2,arg3);
				PRINT_LOG(work);
				if (start_exec(4,arg,0,NULL,&startexec_info.ps,&startexec_info.thread) == -1) {
					retCode = STATUS_EXEC_FAIL;errLine = __LINE__;
					break;
				}
				startexec_info.state = 1;
				WaitForSingleObject(startexec_info.ps,0xFFFFFFFF);
				finish_exec(NULL,&startexec_info.ps,&startexec_info.thread);

				retCode = PLG_InfoAudioData(dec_wav,&inFmt,&inSample,&fileInfo);
				if (retCode != STATUS_SUCCESS) {
					PRINT_LOG("Decode Error");
					retCode = STATUS_FILE_READ_ERR;errLine = __LINE__;
					break;
				}
				if (strcmp(inFmt.fmt,"wav") != 0) {
					retCode = STATUS_FILE_READ_ERR;errLine = __LINE__;
					break;
				}
				// 元の入力ファイル名を保存(後でメタ情報の復元時に使用する)
				fprintf(fp_param,"%s\n",argv[1]);
				strcpy(argv[1],dec_wav);
				fprintf(fp_files,"%s\n",dec_wav);
			} else {
				fprintf(fp_param,"%s\n",argv[1]);
			}

			// r1,r2,r3,r4,r5,r6 ファイルの出力先パスファイル生成
			_splitpath(argv[2],drive,dir,fname,ext);
			_splitpath(workpath,workdrive,workdir,workfname,workext);
			if (strlen(workdrive) == 2 && strlen(workdir) >= 1) {
				strcpy(workfname,fname);
			} else {
				strcpy(workdrive,drive);
				strcpy(workdir,dir);
				strcpy(workfname,fname);
			}
			sprintf(work,"workpath:%s,%s,%s,%s",workdrive,workdir,workfname,workext);
			PRINT_LOG(work);

#ifdef _OPENMP
	omp_set_num_threads(thread);
#endif
			if (strcmp(inFmt.fmt,"dsf") == 0) {
				/* DSF */
				fprintf(stdout,"[dsf2raw]\n");
				fflush(stdout);
				fclose(fp_param);
				fclose(fp_files);
				param_info2.decode = 1;
				param_info2.argv4  = argv[4];
				dsf_decode(argv[1],argv[2],&param_info2);
			} else {
				/* WAV */
				fprintf(stdout,"[wav2raw]\n");
				fflush(stdout);

				// 入力ファイルオープン
				PRINT_LOG(argv[1]);
				fio_open(&fp_r,argv[1],FIO_MODE_R);
				if (fp_r.error) {
					retCode = STATUS_FILE_READ_ERR;errLine = __LINE__;
					break;
				}
				fio_set_memory_limit(&fp_r,20,fio);

				seekPtr = (fio_size)fileInfo.dataOffset;
				fio_seek(&fp_r,seekPtr,SEEK_SET);
				if (fp_r.error) {
					retCode = STATUS_FILE_READ_ERR;errLine = __LINE__;
					break;
				}

				// パラメーターの追加
				sprintf(work," -is:%d -iw:%d -ch:%d",inFmt.sample,inFmt.bitsPerSample,inFmt.channel);
				strcat(argv[4],work);

				max_size = inSample;
				max_size *= sizeof (SSIZE);

				// 出力ファイル名の作成(Left)
				_makepath(tmppath,workdrive,workdir,workfname,"r1");
				fio_open(&fp_w1,tmppath,FIO_MODE_W);
				if (fp_w1.error) {
					retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
					break;
				}
				fprintf(fp_files,"%s\n",tmppath);
				fio_set_maxsize(&fp_w1,max_size);
				fio_set_memory_limit(&fp_w1,20,fio);

				if (inFmt.channel >= 2) {
					// 出力ファイル名の作成(Right)
					_makepath(tmppath,workdrive,workdir,workfname,"r2");
					fio_open(&fp_w2,tmppath,FIO_MODE_W);
					if (fp_w2.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
					fprintf(fp_files,"%s\n",tmppath);
					fio_set_maxsize(&fp_w2,max_size);
					fio_set_memory_limit(&fp_w2,20,fio);
				}
				if (inFmt.channel >= 3) {
					// 出力ファイル名の作成(3)
					_makepath(tmppath,workdrive,workdir,workfname,"r3");
					fio_open(&fp_w3,tmppath,FIO_MODE_W);
					if (fp_w3.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
					fprintf(fp_files,"%s\n",tmppath);
					fio_set_maxsize(&fp_w3,max_size);
					fio_set_memory_limit(&fp_w3,20,fio);
				}
				if (inFmt.channel >= 4) {
					// 出力ファイル名の作成(4)
					_makepath(tmppath,workdrive,workdir,workfname,"r4");
					fio_open(&fp_w4,tmppath,FIO_MODE_W);
					if (fp_w4.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
					fprintf(fp_files,"%s\n",tmppath);
					fio_set_maxsize(&fp_w4,max_size);
					fio_set_memory_limit(&fp_w4,20,fio);
				}
				if (inFmt.channel >= 5) {
					// 出力ファイル名の作成(5)
					_makepath(tmppath,workdrive,workdir,workfname,"r5");
					fio_open(&fp_w5,tmppath,FIO_MODE_W);
					if (fp_w5.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
					fprintf(fp_files,"%s\n",tmppath);
					fio_set_maxsize(&fp_w5,max_size);
					fio_set_memory_limit(&fp_w5,20,fio);
				}
				if (inFmt.channel >= 6) {
					// 出力ファイル名の作成(6)
					_makepath(tmppath,workdrive,workdir,workfname,"r6");
					// 出力ファイルオープン
					fio_open(&fp_w6,tmppath,FIO_MODE_W);
					if (fp_w6.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
					fprintf(fp_files,"%s\n",tmppath);
					fio_set_maxsize(&fp_w6,max_size);
					fio_set_memory_limit(&fp_w6,20,fio);
				}
				
				// 削除ファイル名
				_makepath(tmppath,workdrive,workdir,workfname,"r1.param");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r2.param");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r3.param");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r4.param");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r5.param");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r6.param");
				fprintf(fp_files,"%s\n",tmppath);

				_makepath(tmppath,workdrive,workdir,workfname,"r1.tmp");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r2.tmp");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r3.tmp");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r4.tmp");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r5.tmp");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r6.tmp");
				fprintf(fp_files,"%s\n",tmppath);

				_makepath(tmppath,workdrive,workdir,workfname,"r1.tmp2");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r2.tmp2");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r3.tmp2");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r4.tmp2");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r5.tmp2");
				fprintf(fp_files,"%s\n",tmppath);
				_makepath(tmppath,workdrive,workdir,workfname,"r6.tmp2");
				fprintf(fp_files,"%s\n",tmppath);

				// 読み込みバッファサイズ計算
				if (inFmt.bitsPerSample != 20) {
					rs = (inFmt.bitsPerSample / 8) * inFmt.channel;
				} else {
					rs = (24 / 8) * inFmt.channel;
				}

				inBuffer = (unsigned char *)malloc(rs * inFmt.sample * 10);
				if (inBuffer == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;errLine = __LINE__;
					break;
				}
				// 変換と書き出し

				err1 = err2 = err3 = err4 = err5 = err6 = 0;
				per = -1;

				for (i = 0;i < inSample;i += inFmt.sample * 10) {
					persent = ((double)i / inSample);
					persent *= 100;
					if (persent != per) {
						fprintf(stdout,"%d%%\n",(int)persent);
						fflush(stdout);
					}
					per = persent;

					if (err1 || err2 || err3 || err4 || err5 || err6) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
					nSample = inFmt.sample * 10;
					if (i + nSample > inSample) {
						nSample = inSample - i;
					}
					rd = fio_read(inBuffer,rs,nSample,&fp_r);
					if (fp_r.error) {
						retCode = STATUS_FILE_READ_ERR;errLine = __LINE__;
						break;
					}
					if (rd == 0) {
						break;
					}
					#pragma omp parallel
					{
						#pragma omp sections
						{
							#pragma omp section
							{
								if (nBitTo64S(inFmt.channel,0,inFmt.bitsPerSample,inBuffer,&fp_w1,rd)) {
									err1 = 1;
								}
							}
							#pragma omp section
							{
								if (inFmt.channel >= 2) {
									if (nBitTo64S(inFmt.channel,1,inFmt.bitsPerSample,inBuffer,&fp_w2,rd)) {
										err2 = 1;
									}
								}
							}
							#pragma omp section
							{
								if (inFmt.channel >= 3) {
									if (nBitTo64S(inFmt.channel,2,inFmt.bitsPerSample,inBuffer,&fp_w3,rd)) {
										err3 = 1;
									}
								}
							}
							#pragma omp section
							{
								if (inFmt.channel >= 4) {
									if (nBitTo64S(inFmt.channel,3,inFmt.bitsPerSample,inBuffer,&fp_w4,rd)) {
										err4 = 1;
									}
								}
							}
							#pragma omp section
							{
								if (inFmt.channel >= 5) {
									if (nBitTo64S(inFmt.channel,4,inFmt.bitsPerSample,inBuffer,&fp_w5,rd)) {
										err5 = 1;
									}
								}
							}
							#pragma omp section
							{
								if (inFmt.channel >= 6) {
									if (nBitTo64S(inFmt.channel,5,inFmt.bitsPerSample,inBuffer,&fp_w6,rd)) {
										err6 = 1;
									}
								}
							}
						}
					}
				}
				fio_close(&fp_r);
				fio_close(&fp_w1);
				if (fp_w1.error) {
					retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
					break;
				}
				if (inFmt.channel >= 2) {
					fio_close(&fp_w2);
					if (fp_w2.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
				}
				if (inFmt.channel >= 3) {
					fio_close(&fp_w3);
					if (fp_w3.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
				}
				if (inFmt.channel >= 4) {
					fio_close(&fp_w4);
					if (fp_w4.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
				}
				if (inFmt.channel >= 5) {
					fio_close(&fp_w5);
					if (fp_w5.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
				}
				if (inFmt.channel >= 6) {
					fio_close(&fp_w6);
					if (fp_w6.error) {
						retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
						break;
					}
				}

				if (retCode != STATUS_SUCCESS) {
					break;
				}
				
				// 音の最大レベルを記録する(1)
				if (NormInfo.maxR1 < 0) {
					NormInfo.maxR1 *= -1;
				}
				if (NormInfo.minR1 < 0) {
					NormInfo.minR1 *= -1;
				}
				max = NormInfo.maxR1;
				if (max < NormInfo.minR1) {
					max = NormInfo.minR1;
				}
				if (NormInfo.cntLv1 > 0) {
					NormInfo.tmpLv1 /= NormInfo.cntLv1;
				}

				if (NormInfo.maxLv1 > 0) {
					avg = (NormInfo.maxLv1 + NormInfo.tmpLv1) / 2;
				} else {
					avg = NormInfo.tmpLv1;
				}
				avg <<= 40;
				persent = (double)max / (double)0x7FFFFFFFFFFFFF;
				wr = fprintf(fp_param,"r1=%.10lf,%llx\n",persent,avg);	// 音のレベル(%),音の平均値
				if (wr == EOF) {
					retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
					break;
				}

				// 音の最大レベルを記録する(2)
				persent = 1.00;
				avg = 0;
				if (inFmt.channel >= 2) {
					if (NormInfo.maxR2 < 0) {
						NormInfo.maxR2 *= -1;
					}
					if (NormInfo.minR2 < 0) {
						NormInfo.minR2 *= -1;
					}
					max = NormInfo.maxR2;
					if (max < NormInfo.minR2) {
						max = NormInfo.minR2;
					}

					if (NormInfo.cntLv2 > 0) {
						NormInfo.tmpLv2 /= NormInfo.cntLv2;
					}
					if (NormInfo.maxLv2 > 0) {
						avg = (NormInfo.maxLv2 + NormInfo.tmpLv2) / 2;
					} else {
						avg = NormInfo.tmpLv2;
					}
					avg <<= 40;
					persent = (double)max / (double)0x7FFFFFFFFFFFFF;
				}
				wr = fprintf(fp_param,"r2=%.10lf,%llx\n",persent,avg);
				if (wr == EOF) {
					retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
					break;
				}
				// 音の最大レベルを記録する(3)
				persent = 1.00;
				avg = 0;
				if (inFmt.channel >= 3) {
					if (NormInfo.maxR3 < 0) {
						NormInfo.maxR3 *= -1;
					}
					if (NormInfo.minR3 < 0) {
						NormInfo.minR3 *= -1;
					}
					max = NormInfo.maxR3;
					if (max < NormInfo.minR3) {
						max = NormInfo.minR3;
					}
					if (NormInfo.cntLv3 > 0) {
						NormInfo.tmpLv3 /= NormInfo.cntLv3;
					}
					if (NormInfo.maxLv3 > 0) {
						avg = (NormInfo.maxLv3 + NormInfo.tmpLv3) / 2;
					} else {
						avg = NormInfo.tmpLv3;
					}
					avg <<= 40;
					persent = (double)max / (double)0x7FFFFFFFFFFFFF;
				}
				wr = fprintf(fp_param,"r3=%.10lf,%llx\n",persent,avg);
				if (wr == EOF) {
					retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
					break;
				}
				// 音の最大レベルを記録する(4)
				persent = 1.00;
				avg = 0;
				if (inFmt.channel >= 4) {
					if (NormInfo.maxR4 < 0) {
						NormInfo.maxR4 *= -1;
					}
					if (NormInfo.minR4 < 0) {
						NormInfo.minR4 *= -1;
					}
					max = NormInfo.maxR4;
					if (max < NormInfo.minR4) {
						max = NormInfo.minR4;
					}
					if (NormInfo.cntLv4 > 0) {
						NormInfo.tmpLv4 /= NormInfo.cntLv4;
					}
					if (NormInfo.maxLv4 > 0) {
						avg = (NormInfo.maxLv4 + NormInfo.tmpLv4) / 2;
					} else {
						avg = NormInfo.tmpLv4;
					}
					avg <<= 40;
					persent = (double)max / (double)0x7FFFFFFFFFFFFF;
				}
				wr = fprintf(fp_param,"r4=%.10lf,%llx\n",persent,avg);
				if (wr == EOF) {
					retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
					break;
				}
				// 音の最大レベルを記録する(5)
				persent = 1.00;
				avg = 0;
				if (inFmt.channel >= 5) {
					if (NormInfo.maxR5 < 0) {
						NormInfo.maxR5 *= -1;
					}
					if (NormInfo.minR5 < 0) {
						NormInfo.minR5 *= -1;
					}
					max = NormInfo.maxR5;
					if (max < NormInfo.minR5) {
						max = NormInfo.minR5;
					}
					if (NormInfo.cntLv5 > 0) {
						NormInfo.tmpLv5 /= NormInfo.cntLv5;
					}
					if (NormInfo.maxLv5 > 0) {
						avg = (NormInfo.maxLv5 + NormInfo.tmpLv5) / 2;
					} else {
						avg = NormInfo.tmpLv5;
					}
					avg <<= 40;
					persent = (double)max / (double)0x7FFFFFFFFFFFFF;
				}
				wr = fprintf(fp_param,"r5=%.10lf,%llx\n",persent,avg);
				if (wr == EOF) {
					retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
					break;
				}
				// 音の最大レベルを記録する(6)
				persent = 1.00;
				if (inFmt.channel >= 6) {
					if (NormInfo.maxR6 < 0) {
						NormInfo.maxR6 *= -1;
					}
					if (NormInfo.minR6 < 0) {
						NormInfo.minR6 *= -1;
					}
					max = NormInfo.maxR6;
					if (max < NormInfo.minR6) {
						max = NormInfo.minR6;
					}
					if (NormInfo.cntLv6 > 0) {
						NormInfo.tmpLv6 /= NormInfo.cntLv6;
					}
					if (NormInfo.maxLv6 > 0) {
						avg = (NormInfo.maxLv6 + NormInfo.tmpLv6) / 2;
					} else {
						avg = NormInfo.tmpLv6;
					}
					avg <<= 40;
					persent = (double)max / (double)0x7FFFFFFFFFFFFF;
				}
				wr = fprintf(fp_param,"r6=%.10lf,%llx\n",persent,avg);
				if (wr == EOF) {
					retCode = STATUS_FILE_WRITE_ERR;errLine = __LINE__;
					break;
				}
				fclose(fp_param);
				fclose(fp_files);
			}
			retCode = STATUS_SUCCESS;
		} while (0);
	}
	if (retCode != STATUS_SUCCESS) {
		_splitpath(argv[2],drive,dir,fname,ext);
		_makepath(tmppath,drive,dir,fname,"err");
		ofp = fopen(tmppath,"w");
		if (ofp) {
			switch (retCode) {
				case STATUS_PARAMETER_ERR:
					fprintf(ofp,"to_raw_main:[%04d] Parameter error.\n",errLine);
					break;
				case STATUS_FILE_READ_ERR:
					fprintf(ofp,"to_raw_main:[%04d] File read error.\n",errLine);
					break;
				case STATUS_FILE_WRITE_ERR:
					fprintf(ofp,"to_raw_main:[%04d] File write error.\n",errLine);
					break;
				case STATUS_MEM_ALLOC_ERR:
					fprintf(ofp,"to_raw_main:[%04d] Memory Allocation error.\n",errLine);
					break;
				default:
					fprintf(ofp,"to_raw_main:[%04d] Other error.\n",errLine);
			}
			fclose(ofp);
		}
	}

	return retCode;
}

//---------------------------------------------------------------------------
// Function   : nBitTo64S
// Description: nBit のデータを64Bit 符号付データに変換する(内部表現は60Bit)
// ---
//	nCh		:チャンネル数
//	ch		:音声データのチャネル
//	bit 	:入力データのビット数
//	in		:入力データのアドレス
//	fio 	:出力FIOのアドレス
//	nSample :サンプル数
//
int nBitTo64S(int nCh,int ch,int bit,void *in,FIO *fio,DWORD nSample)
{
	DWORD i;
	SSIZE ns;
	short *p16;
	unsigned char *p24;
	float  *p32;
	double *p64;
	SSIZE out;
	SSIZE max,min;
	SSIZE maxLv,maxLv2;
	int next;
	fio_size ws;

	ns = 0;
	maxLv2 = 0;
	if (bit == 16) {
		/* 16Bit */
		next = 1 * nCh;
		p16 = (short *)in;
		if (ch >= 1) {
			p16++;
		}
		if (ch >= 2) {
			p16++;
		}
		if (ch >= 3) {
			p16++;
		}
		if (ch >= 4) {
			p16++;
		}
		if (ch >= 5) {
			p16++;
		}
	}
	if (bit == 24 || bit == 20) {
		/* 24Bit */
		next = 3 * nCh;
		p24 = (char *)in;
		if (ch >= 1) {
			p24+=3;
		}
		if (ch >= 2) {
			p24+=3;
		}
		if (ch >= 3) {
			p24+=3;
		}
		if (ch >= 4) {
			p24+=3;
		}
		if (ch >= 5) {
			p24+=3;
		}
	}
	if (bit == 32) {
		/* 32Bit */
		p32 = (float *)in;
		next = 1 * nCh;
		if (ch >= 1) {
			p32++;
		}
		if (ch >= 2) {
			p32++;
		}
		if (ch >= 3) {
			p32++;
		}
		if (ch >= 4) {
			p32++;
		}
		if (ch >= 5) {
			p32++;
		}
	}
	if (bit == 64) {
		/* 64Bit */
		p64 = (double *)in;
		next = 1 * nCh;
		if (ch >= 1) {
			p64++;
		}
		if (ch >= 2) {
			p64++;
		}
		if (ch >= 3) {
			p64++;
		}
		if (ch >= 4) {
			p64++;
		}
		if (ch >= 5) {
			p64++;
		}
	}

	if (ch == 0) {
		max = NormInfo.maxR1;
		min = NormInfo.minR1;
		maxLv = NormInfo.maxLv1;
		maxLv2 = NormInfo.tmpLv1;
		ns	   = NormInfo.cntLv1;
	} else if (ch == 1) {
		max = NormInfo.maxR2;
		min = NormInfo.minR2;
		maxLv = NormInfo.maxLv2;
		maxLv2 = NormInfo.tmpLv2;
		ns	   = NormInfo.cntLv2;
	} else if (ch == 2) {
		max = NormInfo.maxR3;
		min = NormInfo.minR3;
		maxLv = NormInfo.maxLv3;
		maxLv2 = NormInfo.tmpLv3;
		ns	   = NormInfo.cntLv3;
	} else if (ch == 3) {
		max = NormInfo.maxR4;
		min = NormInfo.minR4;
		maxLv = NormInfo.maxLv4;
		maxLv2 = NormInfo.tmpLv4;
		ns	   = NormInfo.cntLv4;
	} else if (ch == 4) {
		max = NormInfo.maxR5;
		min = NormInfo.minR5;
		maxLv = NormInfo.maxLv5;
		maxLv2 = NormInfo.tmpLv5;
		ns	   = NormInfo.cntLv5;
	} else if (ch == 5) {
		max = NormInfo.maxR6;
		min = NormInfo.minR6;
		maxLv = NormInfo.maxLv6;
		maxLv2 = NormInfo.tmpLv6;
		ns	   = NormInfo.cntLv6;
	}

	switch (bit) {
		case 16:
			for (i = 0;i < nSample;i++) {
				out = (SSIZE)*p16;
				out <<= (64 - 16);
				out >>= 8;
				
				if (max < out) {
					max = out;
				} else if (min > out) {
					min = out;
				}
				if ((out >> 40) > 0) {
					maxLv2 += (out >> 40);
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
				p16 += next;
				ws = fio_write(&out,sizeof (SSIZE),1,fio);
				if (fio->error || ws != 1) {
					break;
				}
			}
			break;
		case 20:
			for (i = 0;i < nSample;i++) {
				out  = (SSIZE)p24[2];
				out  = out << 8;
				out |= (SSIZE)p24[1];
				out  = out << 8;
				out |= (SSIZE)p24[0];
				out  = out << (64 - 20);
				out  = out >> 8;

				if (max < out) {
					max = out;
				} else if (min > out) {
					min = out;
				}
				if ((out >> 40) > 0) {
					maxLv2 += (out >> 40);
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

				p24 += next;
				ws = fio_write(&out,sizeof (SSIZE),1,fio);
				if (fio->error || ws != 1) {
					break;
				}
			}
			break;
		case 24:
			for (i = 0;i < nSample;i++) {
				out  = (SSIZE)p24[2];
				out  = out << 8;
				out |= (SSIZE)p24[1];
				out  = out << 8;
				out |= (SSIZE)p24[0];
				out  = out << (64 - 24);
				out  = out >> 8;

				if (max < out) {
					max = out;
				} else if (min > out) {
					min = out;
				}
				if ((out >> 40) > 0) {
					maxLv2 += (out >> 40);
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

				p24 += next;
				ws = fio_write(&out,sizeof (SSIZE),1,fio);
				if (fio->error || ws != 1) {
					break;
				}
			}
			break;
		case 32:
			for (i = 0;i < nSample;i++) {
				out = (SSIZE)(*p32 * 0x7FFFFFFFFFFFFF);
				
				if (max < out) {
					max = out;
				} else if (min > out) {
					min = out;
				}
				if ((out >> 40) > 0) {
					maxLv2 += (out >> 40);
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
				
				p32 += next;
				ws = fio_write(&out,sizeof (SSIZE),1,fio);
				if (fio->error || ws != 1) {
					break;
				}
			}
			break;
		case 64:
			for (i = 0;i < nSample;i++) {
				out = (SSIZE)(*p64 * 0x7FFFFFFFFFFFFF);
				
				if (max < out) {
					max = out;
				} else if (min > out) {
					min = out;
				}
				if ((out >> 40) > 0) {
					maxLv2 += (out >> 40);
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

				p64 += next;
				ws = fio_write(&out,sizeof (SSIZE),1,fio);
				if (fio->error || ws != 1) {
					break;
				}
			}
			break;
		default:
			break;
	}
	
	if (ch == 0) {
		NormInfo.maxR1 = max;
		NormInfo.minR1 = min;
		NormInfo.maxLv1 = maxLv;
		NormInfo.tmpLv1 = maxLv2;
		NormInfo.cntLv1 = ns;
	} else if (ch == 1) {
		NormInfo.maxR2 = max;
		NormInfo.minR2 = min;
		NormInfo.maxLv2 = maxLv;
		NormInfo.tmpLv2 = maxLv2;
		NormInfo.cntLv2 = ns;
	} else if (ch == 2) {
		NormInfo.maxR3 = max;
		NormInfo.minR3 = min;
		NormInfo.maxLv3 = maxLv;
		NormInfo.tmpLv3 = maxLv2;
		NormInfo.cntLv3 = ns;
	} else if (ch == 3) {
		NormInfo.maxR4 = max;
		NormInfo.minR4 = min;
		NormInfo.maxLv4 = maxLv;
		NormInfo.tmpLv4 = maxLv2;
		NormInfo.cntLv4 = ns;
	} else if (ch == 4) {
		NormInfo.maxR5 = max;
		NormInfo.minR5 = min;
		NormInfo.maxLv5 = maxLv;
		NormInfo.tmpLv5 = maxLv2;
		NormInfo.cntLv5 = ns;
	} else if (ch == 5) {
		NormInfo.maxR6 = max;
		NormInfo.minR6 = min;
		NormInfo.maxLv6 = maxLv;
		NormInfo.tmpLv6 = maxLv2;
		NormInfo.cntLv6 = ns;
	}
	if (fio->error || ws != 1) {
		return -1;
	}
	
	return 0;
}

//---------------------------------------------------------------------------
// Function   : to_wav_main
// Description: wav ファイル化メイン
//
//
int to_wav_main(int argc, char *argv[])
{
	char workpath[_MAX_PATH];
	char tmppath[_MAX_PATH];
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char tmpfn[_MAX_PATH];
	char ext[_MAX_EXT];
	static char workdrive[_MAX_DRIVE];
	static char workdir[_MAX_DIR];
	static char workfname[_MAX_FNAME];
	static char workext[_MAX_EXT];
	static char pdrive[_MAX_DRIVE];
	static char pdir[_MAX_DIR];
	static char pfname[_MAX_FNAME];
	static char pext[_MAX_EXT];
	char *arg_encorder[3];
	char work[2048];
	char param[2048];
	char fn[6][10]={"_1","_2","_3","_4","_5","_6"};
	char s[300];
	int paramFlag;
	int ccc;
	FILE *fp;
	FIO fp_r1,fp_r2,fp_r3,fp_r4,fp_r5,fp_r6,fp_w;
	FIO fp_ws[6];
	FIO *p_fp1,*p_fp2,*p_fp3,*p_fp4,*p_fp5,*p_fp6;
	SOUNDFMT inFmt;
	SOUNDFMT outFmt;
	SOUNDFMT flacFmt;
	char *p1,*p2;
	DWORD inSample,outSample,flacSample,startSample;
	SSIZE size;
	double nx;
	long rd;
	int i,count;
	int retCode;
	int thread;
	int fio;
	int genCh;
	int orgfile_is_dsf;
	double *p_nx[6];
	SSIZE  *p_max[6];
	SSIZE  *p_avg[6];
	FILEINFO fileInfo;
	long temp,temp2,temp3,temp4;
	double perR1,perR2,perR3,perR4,perR5,perR6;
	SSIZE maxR1,maxR2,maxR3,maxR4,maxR5,maxR6;
	SSIZE avgR1,avgR2,avgR3,avgR4,avgR5,avgR6;
	double per;
	SSIZE max,avg;
	SSIZE max_l,max_r;
	thread = 1;
	genCh = 0;
	orgfile_is_dsf = 0;
	do {
		param[0] = '\0';
		memset(&paramInfo,0,sizeof (PARAM_INFO));
		memset(&fileInfo,0,sizeof (FILEINFO));
		memset(&fp_ws,0,sizeof (FIO) * 6);
		paramFlag = 0;
		paramInfo.chC = 0;
		paramInfo.chS = 0;
		paramInfo.chLFE = 0;
		paramInfo.fio = 5;
		paramInfo.dsd_fmt = -1;
		paramInfo.mode = 0;
		paramInfo.norm = 0;
		paramInfo.norm_option = 1;
		paramInfo.ditherLv = 0;
		paramInfo.encorder = 0;
		arg_encorder[0] = NULL;
		arg_encorder[1] = NULL;
		arg_encorder[2] = NULL;
		// Wave ファイルのは違う順番
		p_nx[0]  = &paramInfo.nx1;	// Left
		p_nx[1]  = &paramInfo.nx2;	// Right
		p_nx[2]  = &paramInfo.nx3;	// Center
		p_nx[3]  = &paramInfo.nx4;	// Surround Left
		p_nx[4]  = &paramInfo.nx5;	// Surround Right
		p_nx[5]  = &paramInfo.nx6;	// LFE
		p_max[0] = &maxR1;
		p_max[1] = &maxR2;
		p_max[2] = &maxR3;
		p_max[3] = &maxR4;
		p_max[4] = &maxR5;
		p_max[5] = &maxR6;
		p_avg[0] = &avgR1;
		p_avg[1] = &avgR2;
		p_avg[2] = &avgR3;
		p_avg[3] = &avgR4;
		p_avg[4] = &avgR5;
		p_avg[5] = &avgR6;

		p_fp1 = p_fp2 = p_fp3 = p_fp4 = p_fp5 = p_fp6 = NULL;
		
		retCode = STATUS_SUCCESS;
		paramInfo.errLine = __LINE__;
		PRINT_LOG("before argc=5");

		if (argc == 5) {
//			strcpy(paramInfo.fromfile,argv[1]);

			// default parameter
			fp = fopen(argv[3],"r");
			if (fp == NULL) {
				retCode = STATUS_PARAMETER_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fgets(work,2047,fp) == NULL) {
				retCode = STATUS_PARAMETER_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			p1 = strrchr(work,'\n');if (p1 != NULL) *p1 = '\0';
			strcat(param,work);
			strcat(param," ");
			if (strlen(argv[4]) >= 1) strcat(param,argv[4]);

			// テンポラリファイルの出力先
			if (fgets(workpath,_MAX_PATH - 1,fp) == NULL) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			p1 = strrchr(workpath,'\n');if (p1 != NULL) *p1 = '\0';
			if (strlen(workpath) >= 2 && workpath[strlen(workpath) - 1] != '\\') strcat(workpath,"\\");

			// FLAC(Encoder Option)
			if (fgets(work,500,fp) == NULL) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			p1 = strrchr(work,'\n');if (p1 != NULL) *p1 = '\0';
			strcpy(paramInfo.opt_flac,"|");
			strcat(paramInfo.opt_flac,work);
			
			// FLAC(WavPack Option)
			if (fgets(work,500,fp) == NULL) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			p1 = strrchr(work,'\n');if (p1 != NULL) *p1 = '\0';
			strcpy(paramInfo.opt_wavpack,"|");
			strcat(paramInfo.opt_wavpack,work);

			// FLAC(WavPack Option)
			if (fgets(work,500,fp) == NULL) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			p1 = strrchr(work,'\n');if (p1 != NULL) *p1 = '\0';
			strcpy(paramInfo.opt_mp3,"|");
			strcat(paramInfo.opt_mp3,work);

			// param ファイル
			PRINT_LOG("before open param");
			_splitpath(argv[2],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,fname,"param");
			// ファイルオープン
			fp = fopen(tmppath,"r");
			if (fp == NULL) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			
			if (fgets(work,_MAX_PATH,fp) == NULL) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			p1 = strrchr(work,'\n');if (p1 != NULL) *p1 = '\0';
			retCode = PLG_InfoAudioData(work,&inFmt,&inSample,&fileInfo);
			if (retCode == STATUS_SUCCESS) {
				if (strcmp(inFmt.fmt,"dsf") == 0) {
					orgfile_is_dsf = 1;
				}
			}
			strcpy(paramInfo.fromfile,work);

			perR1 = perR2 = perR3 = perR4 = perR5 = perR6 = 0;
			PRINT_LOG("before r1");
			if (fscanf(fp,"r1=%lf,%llx\n",&perR1,&avgR1) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			PRINT_LOG("before r2");
			if (fscanf(fp,"r2=%lf,%llx\n",&perR2,&avgR2) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fscanf(fp,"r3=%lf,%llx\n",&perR3,&avgR3) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fscanf(fp,"r4=%lf,%llx\n",&perR4,&avgR4) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fscanf(fp,"r5=%lf,%llx\n",&perR5,&avgR5) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fscanf(fp,"r6=%lf,%llx\n",&perR6,&avgR6) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			PRINT_LOG("r1 r2");

			maxR1 = (SSIZE)(perR1 * 0x7FFFFFFFFFFFFF);
			maxR2 = (SSIZE)(perR2 * 0x7FFFFFFFFFFFFF);
			maxR3 = (SSIZE)(perR3 * 0x7FFFFFFFFFFFFF);
			maxR4 = (SSIZE)(perR4 * 0x7FFFFFFFFFFFFF);
			maxR5 = (SSIZE)(perR5 * 0x7FFFFFFFFFFFFF);
			maxR6 = (SSIZE)(perR6 * 0x7FFFFFFFFFFFFF);

			fclose(fp);
			p1 = param;
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
							paramInfo.sampling = temp;
							paramFlag |= 0x01;
							break;
					}
				}
				if (sscanf(p1,"-w:%ld",&temp) == 1) {
					switch (temp) {
						case 16:
						case 24:
						case 32:
						case 64:
							paramInfo.bitwidth = temp;
							paramFlag |= 0x02;
							break;
					}
				}
				if (sscanf(p1,"-ch:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 6) {
						paramInfo.ch = temp;
						paramFlag |= 0x04;
					}
				}

				if (sscanf(p1,"-output_dither:%ld",&temp) == 1) {
					switch (temp) {
						case 0:
						case 1:
						case 2:
						case 3:
							paramInfo.mode = temp;
							break;
					}
				}
				if (sscanf(p1,"-output_dither_option:%ld",&temp) == 1) {
					if (temp >= 0 && temp <= 16) {
						paramInfo.ditherLv = temp;
					}
				}
				if (sscanf(p1,"-norm:%ld",&temp) == 1) {
					switch (temp) {
						case 0:
						case 1:
							paramInfo.norm = temp;
							break;
					}
				}
				if (sscanf(p1,"-norm_option:%ld",&temp) == 1) {
					switch (temp) {
						case 1:
						case 2:
							paramInfo.norm_option = temp;
							break;
					}
				}
				if (strcmp(p1,"-bwf") == 0) {
					paramInfo.bwf = 1;
				}
				if (strcmp(p1,"-raw") == 0) {
					paramInfo.raw = 1;
				}
				if (strcmp(p1,"-rf64") == 0) {
					paramInfo.rf64 = 1;
				}
				if (strcmp(p1,"-C") == 0) {
					paramInfo.chC = 1;
					genCh++;
				}
				if (strcmp(p1,"-SLR") == 0) {
					paramInfo.chS = 1;
					genCh += 2;
				}
				if (strcmp(p1,"-LFE") == 0) {
					paramInfo.chLFE = 1;
					genCh++;
				}
				if (sscanf(p1,"-MC_Option:%d,%d,%d,%d",&temp,&temp2,&temp3,&temp4) == 4) {
					if (temp == 1) {
						paramInfo.split = 1;
					}
					if (temp2 == 1) {
						paramInfo.downmix = 1;
					}
				}
				if (sscanf(p1,"-hfa:%ld",&temp) == 1) {
					switch (temp) {
						case 0:
						case 1:
						case 2:
						case 3:
							paramInfo.hfa = temp;
							break;
					}
				}
				if (strcmp(p1,"-abe") == 0) {
					paramInfo.abe = 1;
				}
				if (sscanf(p1,"-encorder:%ld",&temp) == 1) {
					if (temp >= 0 && temp <= 3) paramInfo.encorder = temp;
				}
				if (sscanf(p1,"-thread:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 24) {
						thread = (int)temp;
					}
				}
				if (sscanf(p1,"-fio:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 16000) {
						paramInfo.fio = temp;
					}
				}
				if (sscanf(p1,"-dsd_fmt:%ld",&temp) == 1) {
					if (temp == 64 || temp == 128 || temp == 256) {
						paramInfo.dsd_fmt = (int)temp;
					}
				}

				if (p2 == NULL) {
					break;
				}
				p1 = p2 + 1;
				p2 = strchr(p1,(int)' ');
			}
			if (orgfile_is_dsf == 1) {
				paramInfo.norm = 1;
				paramInfo.norm_option = 1;
			}
			if (paramFlag == 0x07) {
#ifdef _OPENMP
	omp_set_num_threads(thread);
#endif
				if (paramInfo.dsd_fmt != -1) {
					paramInfo.downmix = 1;
					paramInfo.split   = 0;
				}

				if (paramInfo.downmix) {
					genCh = 0;
				}
				paramInfo.ch += genCh;

				if (paramInfo.encorder == 1) {
					_splitpath(argv[0],pdrive,pdir,pfname,pext);
					_makepath(tmppath,pdrive,pdir,"flac","exe");
					arg_encorder[0] = strdup(tmppath);
					arg_encorder[1] = strdup(paramInfo.opt_flac);
					arg_encorder[2] = malloc(4096);
					if (arg_encorder[0] == NULL || arg_encorder[1] == NULL || arg_encorder[2] == NULL) {
						retCode = STATUS_MEM_ALLOC_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					strcpy(arg_encorder[2],"");
				} else if (paramInfo.encorder == 2) {
					_splitpath(argv[0],pdrive,pdir,pfname,pext);
					_makepath(tmppath,pdrive,pdir,"wavpack","exe");
					arg_encorder[0] = strdup(tmppath);
					arg_encorder[1] = strdup(paramInfo.opt_wavpack);
					arg_encorder[2] = malloc(4096);
					if (arg_encorder[0] == NULL || arg_encorder[1] == NULL || arg_encorder[2] == NULL) {
						retCode = STATUS_MEM_ALLOC_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					strcpy(arg_encorder[2],"");
				} else if (paramInfo.encorder == 3) {
					_splitpath(argv[0],pdrive,pdir,pfname,pext);
					_makepath(tmppath,pdrive,pdir,"lame","exe");
					arg_encorder[0] = strdup(tmppath);
					arg_encorder[1] = strdup(paramInfo.opt_mp3);
					arg_encorder[2] = malloc(4096);
					if (arg_encorder[0] == NULL || arg_encorder[1] == NULL || arg_encorder[2] == NULL) {
						retCode = STATUS_MEM_ALLOC_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					strcpy(arg_encorder[2],"");
				}

				retCode = PLG_InfoAudioData(argv[1],&inFmt,&inSample,&fileInfo);
				if (retCode != STATUS_SUCCESS) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}

				outFmt.sample  = paramInfo.sampling;
				if (paramInfo.split == 0) {
					outFmt.channel = paramInfo.ch;
				} else {
					outFmt.channel = 1;
				}
				outFmt.bitsPerSample = (unsigned char)paramInfo.bitwidth;

				if (paramInfo.ch == 3) {
					if (paramInfo.chC == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_max[2] = &maxR3;
						p_avg[2] = &avgR3;
						strcpy(fn[2],"_C");
					} else if (paramInfo.chLFE == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_max[2] = &maxR6;
						p_avg[2] = &avgR6;
						strcpy(fn[2],"_LFE");
					}
				}
				if (paramInfo.ch == 4) {
					if (paramInfo.chS == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_nx[3]  = &paramInfo.nx4;
						p_max[2] = &maxR4;
						p_max[3] = &maxR5;
						p_avg[2] = &avgR4;
						p_avg[3] = &avgR5;
						strcpy(fn[2],"_SL");
						strcpy(fn[3],"_SR");
					}
					if (paramInfo.chC == 1 && paramInfo.chLFE == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_nx[3]  = &paramInfo.nx4;
						p_max[2] = &maxR3;
						p_max[3] = &maxR6;
						p_avg[2] = &avgR3;
						p_avg[3] = &avgR6;
						strcpy(fn[2],"_C");
						strcpy(fn[3],"_LFE");
					}
				}
				if (paramInfo.ch == 5) {
					if (paramInfo.chC == 1 && paramInfo.chS == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_nx[3]  = &paramInfo.nx4;
						p_nx[4]  = &paramInfo.nx5;
						p_max[2] = &maxR3;
						p_max[3] = &maxR4;
						p_max[4] = &maxR5;
						p_avg[2] = &avgR3;
						p_avg[3] = &avgR4;
						p_avg[4] = &avgR5;
						strcpy(fn[2],"_C");
						strcpy(fn[3],"_SL");
						strcpy(fn[4],"_SR");
					}
					if (paramInfo.chLFE == 1 && paramInfo.chS == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_nx[3]  = &paramInfo.nx4;
						p_nx[4]  = &paramInfo.nx5;
						p_max[2] = &maxR6;
						p_max[3] = &maxR4;
						p_max[4] = &maxR5;
						p_avg[2] = &avgR6;
						p_avg[3] = &avgR4;
						p_avg[4] = &avgR5;
						strcpy(fn[2],"_LFE");
						strcpy(fn[3],"_SL");
						strcpy(fn[4],"_SR");
					}
				}
				if (paramInfo.ch == 6) {
					p_nx[2]  = &paramInfo.nx3;
					p_nx[5]  = &paramInfo.nx6;
					p_nx[3]  = &paramInfo.nx4;
					p_nx[4]  = &paramInfo.nx5;
					strcpy(fn[2],"_C");
					strcpy(fn[3],"_LFE");
					strcpy(fn[4],"_SL");
					strcpy(fn[5],"_SR");
				}

				// r1 r2 r3 r4 r5 r6ファイルのパス生成
				_splitpath(argv[2],drive,dir,fname,ext);
				_splitpath(workpath,workdrive,workdir,workfname,workext);
				if (strlen(workdrive) == 2 && strlen(workdir) >= 1) {
					strcpy(workfname,fname);
					strcpy(workext,ext);
				} else {
					strcpy(workdrive,drive);
					strcpy(workdir,dir);
					strcpy(workfname,fname);
				}

				_makepath(tmppath,workdrive,workdir,workfname,"r1.param");
				fp = fopen(tmppath,"rb");
				if (fp == NULL) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}
				rd = fscanf(fp,"%lf,%llx",&per,&avg);
				if (rd != 2) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}
				fclose(fp);
				nx = 0;
				if (avg > 10000000000) {
					nx = (double)*p_avg[0] / avg;
					if (nx > 0) {
						max = (SSIZE)((double)*p_max[0] / nx);
					}
				} else {
					max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
					max = (SSIZE)((double)max * 1.02);
				}
				max_l = max;
				nx = 0;
				if (max > 0) {
					if (paramInfo.norm == 1) {
						nx = (double)0x7FFFFFFFFFFFFF / max;
					} else {
						nx = (double)*p_max[0] / max;
					}
					nx *= 0.98;
				}
				*p_nx[0] = nx;

				if (paramInfo.ch >= 2) {
					_makepath(tmppath,workdrive,workdir,workfname,"r2.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[1] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[1] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					max_r = max;
					if (paramInfo.norm == 1) {
						if (paramInfo.norm_option == 1) {
							if (max_l > max_r) {
								max = max_l;
							}
						}
					}
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[1] / max;
						}
						nx *= 0.98;
					}
					*p_nx[1] = nx;
				}

				if (paramInfo.ch == 3) {
					if (paramInfo.chC == 1) {
						_makepath(tmppath,workdrive,workdir,workfname,"r3.param");
					} else if (paramInfo.chLFE) {
						_makepath(tmppath,workdrive,workdir,workfname,"r6.param");
					} else {	// not genCh
						_makepath(tmppath,workdrive,workdir,workfname,"r3.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;
				}
				if (paramInfo.ch == 4 && paramInfo.chS == 1) {
					_makepath(tmppath,workdrive,workdir,workfname,"r4.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					_makepath(tmppath,workdrive,workdir,workfname,"r5.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;
				} else if (paramInfo.ch == 4 && paramInfo.chC == 1 && paramInfo.chLFE == 1) {
					_makepath(tmppath,workdrive,workdir,workfname,"r3.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					_makepath(tmppath,workdrive,workdir,workfname,"r6.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;
				} else if (paramInfo.ch == 4) {
					_makepath(tmppath,workdrive,workdir,workfname,"r3.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					_makepath(tmppath,workdrive,workdir,workfname,"r4.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;
				}

				if (paramInfo.ch == 5) {
					if (paramInfo.chC == 1) {
						_makepath(tmppath,workdrive,workdir,workfname,"r3.param");
					} else if (paramInfo.chLFE) {
						_makepath(tmppath,workdrive,workdir,workfname,"r6.param");
					} else {	// not genCh
						_makepath(tmppath,workdrive,workdir,workfname,"r3.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					_makepath(tmppath,workdrive,workdir,workfname,"r4.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;
					_makepath(tmppath,workdrive,workdir,workfname,"r5.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[4] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[4] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[4] / max;
						}
						nx *= 0.98;
					}
					*p_nx[4] = nx;
				}

				if (paramInfo.ch == 6) {
					_makepath(tmppath,workdrive,workdir,workfname,"r3.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					if (genCh > 0) {
						_makepath(tmppath,workdrive,workdir,workfname,"r6.param");
					} else {
						_makepath(tmppath,workdrive,workdir,workfname,"r4.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;

					if (genCh > 0) {
						_makepath(tmppath,workdrive,workdir,workfname,"r4.param");
					} else {
						_makepath(tmppath,workdrive,workdir,workfname,"r5.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[4] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[4] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[4] / max;
						}
						nx *= 0.98;
					}
					*p_nx[4] = nx;

					if (genCh > 0) {
						_makepath(tmppath,workdrive,workdir,workfname,"r5.param");
					} else {
						_makepath(tmppath,workdrive,workdir,workfname,"r6.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[5] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[5] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[5] / max;
						}
						nx *= 0.98;
					}
					*p_nx[5] = nx;
				}

				_makepath(tmppath,workdrive,workdir,workfname,"r1");
				// ファイルオープン
				fio_open(&fp_r1,tmppath,FIO_MODE_R);
				if (fp_r1.error) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}
				p_fp1 = &fp_r1;
				if (paramInfo.ch >= 2) {
					// R2 ファイル
					_makepath(tmppath,workdrive,workdir,workfname,"r2");

					// ファイルオープン
					fio_open(&fp_r2,tmppath,FIO_MODE_R);
					if (fp_r2.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp2 = &fp_r2;
				}
				if (paramInfo.ch >= 3 && (paramInfo.chC == 1 || genCh == 0)) {
					// R3 ファイル
					_makepath(tmppath,workdrive,workdir,workfname,"r3");

					// ファイルオープン
					fio_open(&fp_r3,tmppath,FIO_MODE_R);
					if (fp_r3.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp3 = &fp_r3;
				}

				if (paramInfo.ch >= 4 && (paramInfo.chS == 1 || genCh == 0)) {
					// R4 ファイル
					_makepath(tmppath,workdrive,workdir,workfname,"r4");

					// ファイルオープン
					fio_open(&fp_r4,tmppath,FIO_MODE_R);
					if (fp_r4.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp4 = &fp_r4;

					// R5 ファイル
					_makepath(tmppath,workdrive,workdir,workfname,"r5");

					// ファイルオープン
					fio_open(&fp_r5,tmppath,FIO_MODE_R);
					if (fp_r5.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp5 = &fp_r5;

				}

				if ((paramInfo.ch >= 3 && paramInfo.chLFE == 1) || (paramInfo.ch == 6 && genCh == 0)) {
					// R6 ファイル
					_makepath(tmppath,workdrive,workdir,workfname,"r6");

					// ファイルオープン
					fio_open(&fp_r6,tmppath,FIO_MODE_R);
					if (fp_r6.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp6 = &fp_r6;
				}

				outSample = 0;
				fio_get_filesize(&fp_r1,&size);
				if (fp_r1.error || size < sizeof (SSIZE)) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}
				outSample = size / sizeof (SSIZE);

				fprintf(stdout,"[raw2wav]\n");
				fflush(stdout);
				if (paramInfo.split == 0) {
					// 出力ファイル
					_splitpath(argv[2],drive,dir,fname,ext);
					if (paramInfo.dsd_fmt == -1) {
						if (paramInfo.raw) {
							_makepath(tmppath,drive,dir,fname,"raw");
						} else {
							_makepath(tmppath,drive,dir,fname,"wav");
						}
						strcpy(paramInfo.tofile,fname);
						strcat(paramInfo.tofile,".wav");
					} else {
						if (paramInfo.dsd_fmt == 64) {
							outFmt.sample = 2822400;
						} else if (paramInfo.dsd_fmt == 128) {
							outFmt.sample = 2822400 * 2;
						} else if (paramInfo.dsd_fmt == 256) {
							outFmt.sample = 2822400 * 4;
						}
						_makepath(tmppath,drive,dir,fname,"dsf");
						strcpy(paramInfo.tofile,fname);
						strcat(paramInfo.tofile,".dsf");
					}
					// ファイルオープン
					fio_open(&fp_w,tmppath,FIO_MODE_W);
					if (fp_w.error) {
						retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fio_set_memory_limit(&fp_w,20,paramInfo.fio);
					paramInfo.fp_w[0] = &fp_w;
					for (count = 1,startSample = 0;count > 0;) {
						if (paramInfo.dsd_fmt == -1) {
							retCode = Normalize(&count,&inFmt,&outFmt,p_fp1,p_fp2,p_fp3,p_fp4,p_fp5,p_fp6,&startSample,outSample,&paramInfo);
						} else {
							PRINT_LOG("");
							retCode = Normalize_DSD(&inFmt,&outFmt,p_fp1,p_fp2,outSample,&paramInfo);
							count = 0;
						}
						fio_close(paramInfo.fp_w[0]);
						if (paramInfo.fp_w[0]->error) {
							break;
						}
						if (paramInfo.encorder != 0) {
							HANDLE process_id,thread_id;
							_splitpath(argv[2],drive,dir,fname,ext);
							_makepath(arg_encorder[2],drive,dir,paramInfo.tofile,"");
							PRINT_LOG("Encorder:");
							PRINT_LOG(arg_encorder[0]);
							PRINT_LOG(arg_encorder[1]);
							PRINT_LOG(arg_encorder[2]);
							if (start_exec(3,arg_encorder,0,NULL,&process_id,&thread_id) == 0) {
								WaitForSingleObject(process_id,0xFFFFFFFF);
								CloseHandle(process_id);
								CloseHandle(thread_id);
							}
						}
						if (count == 0) {
							break;
						}
						count++;
						// 出力ファイル
						_splitpath(argv[2],drive,dir,fname,ext);
						sprintf(fname,"%s_%d.wav",fname,count);
						_makepath(tmppath,drive,dir,fname,NULL);
						// ファイルオープン
						memset(&fp_w,0,sizeof (FIO));
						fio_open(&fp_w,tmppath,FIO_MODE_W);
						if (fp_w.error) {
							retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
							break;
						}
						fio_set_memory_limit(&fp_w,20,paramInfo.fio);
						strcpy(paramInfo.tofile,fname);
						paramInfo.fp_w[0] = &fp_w;
					}
				} else {
					// 出力ファイル
					for (ccc = 0;ccc < paramInfo.ch;ccc++) {
						_splitpath(argv[2],drive,dir,fname,ext);
						strcat(fname,fn[ccc]);
						if (paramInfo.raw) {
							_makepath(tmppath,drive,dir,fname,"raw");
						} else {
							_makepath(tmppath,drive,dir,fname,"wav");
						}
						if (ccc == 0) {
							strcpy(paramInfo.tofile,fname);
							strcat(paramInfo.tofile,"wav");
						}
						// ファイルオープン
						fio_open(&fp_ws[ccc],tmppath,FIO_MODE_W);
						if (fp_ws[ccc].error) {
							retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
							break;
						}
						fio_set_memory_limit(&fp_ws[ccc],20,paramInfo.fio);
						paramInfo.fp_w[ccc] = &fp_ws[ccc];
					}
					for (count = 1,startSample = 0;count > 0;) {
						retCode = Normalize(&count,&inFmt,&outFmt,p_fp1,p_fp2,p_fp3,p_fp4,p_fp5,p_fp6,&startSample,outSample,&paramInfo);
						if (retCode != STATUS_SUCCESS) {
							break;
						}
						for (ccc = 0;ccc < paramInfo.ch;ccc++) {
							fio_close(paramInfo.fp_w[ccc]);
							if (paramInfo.fp_w[ccc]->error) {
								retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
								count = 0;
								break;
							}
						}
						if (count == 0) {
							break;
						}
						count++;
						for (ccc = 0;ccc < paramInfo.ch;ccc++) {
							// 出力ファイル
							_splitpath(argv[2],drive,dir,fname,ext);
							sprintf(fname,"%s_%d%s",fname,count,fn[ccc]);
							_makepath(tmppath,drive,dir,fname,"wav");
							// ファイルオープン
							fio_open(&fp_w,tmppath,FIO_MODE_W);
							if (fp_w.error) {
								retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
								break;
							}
							fio_set_memory_limit(&fp_w,20,paramInfo.fio);
							strcpy(paramInfo.tofile,fname);
							paramInfo.fp_w[ccc] = &fp_w;
							memset(&fp_w,0,sizeof (FIO));
						}
					}
				}
				fio_close(&fp_r1);
				if (p_fp2 != NULL) {
					fio_close(&fp_r2);
				}
				if (p_fp3 != NULL) {
					fio_close(&fp_r3);
				}
				if (p_fp4 != NULL) {
					fio_close(&fp_r4);
				}
				if (p_fp5 != NULL) {
					fio_close(&fp_r5);
				}
				if (p_fp6 != NULL) {
					fio_close(&fp_r6);
				}
			}
		}
		if (!(argc == 5) || paramFlag != 0x07) {
			exit(0);
		}
	} while (0);

	if (arg_encorder[0] != NULL) free(arg_encorder[0]);
	if (arg_encorder[1] != NULL) free(arg_encorder[1]);
	if (arg_encorder[2] != NULL) free(arg_encorder[2]);

	if (retCode != STATUS_SUCCESS) {
		_splitpath(argv[2],drive,dir,fname,ext);
		_makepath(tmppath,drive,dir,fname,"err");
		fp = fopen(tmppath,"a");
		if (fp) {
			switch (retCode) {
				case STATUS_PARAMETER_ERR:
					fprintf(fp,"to_wav_main:[%04d] File read error.\n",paramInfo.errLine);
					break;
				case STATUS_FILE_READ_ERR:
					fprintf(fp,"to_wav_main:[%04d] File read error.\n",paramInfo.errLine);
					break;
				case STATUS_FILE_WRITE_ERR:
					fprintf(fp,"to_wav_main:[%04d] File write error.\n",paramInfo.errLine);
					break;
				case STATUS_MEM_ALLOC_ERR:
					fprintf(fp,"to_wav_main:[%04d] Memory Allocation error.\n",paramInfo.errLine);
					break;
				default:
					fprintf(fp,"to_wav_main:[%04d] Other error.\n",paramInfo.errLine);
			}
			fclose(fp);
		}
	}

	return retCode;
}

//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理
// ---
//	pcount	:ファイル番号のアドレス
//	inFmt	:入力ファイル音声形式情報
//	outFmt	:出力ファイル音声形式情報
//	fp_r1	:音声データのFIO構造体
//	fp_r2	:音声データのFIO構造体
//	fp_r3	:音声データのFIO構造体
//	fp_r4	:音声データのFIO構造体
//	fp_r5	:音声データのFIO構造体
//	fp_r6	:音声データのFIO構造体
//	startSample : 開始サンプル
//	nSample :処理をするサンプル数のアドレス
//	param	:パラメーター構造体
//
int Normalize(int *pCount,SOUNDFMT *inFmt,SOUNDFMT *outFmt,FIO *fp_r1,FIO *fp_r2,FIO *fp_r3,FIO *fp_r4,FIO *fp_r5,FIO *fp_r6,DWORD *startSample,DWORD nSample,PARAM_INFO *param)
{
	DWORD i,j,outSample;
	DWORD chMask;
	BYTE *header;
	BYTE *bwf_chunk;
	BYTE *buf1,*buf2,*buf3,*buf4,*buf5,*buf6;
	FIO *p_fp1,*p_fp2,*p_fp3,*p_fp4,*p_fp5,*p_fp6;
	long header_size;
	fio_size write_size,ws,data_size,file_size;
	fio_size data_start,data_end;
	BROADCAST_EXT *bext;
	BYTE *chunk_data;
	long  chunk_size;
	int retCode,retCode1,retCode2,retCode3,retCode4,retCode5,retCode6;
	int nChunk;
	int bwf_enable;
	long bwf_size;
	int div_flag;
	double persent,per;
	char ppp[50];
	int ccc;
	
	do {
		div_flag = 0;
		header = NULL;
		buf1 = NULL;
		buf2 = NULL;
		buf3 = NULL;
		buf4 = NULL;
		buf5 = NULL;
		buf6 = NULL;
		p_fp1 = p_fp2 = p_fp3 = p_fp4 = p_fp5 = p_fp6 = NULL;

		retCode1 = retCode2 = retCode3 = retCode4 = retCode5 = retCode6 = 0;
		retCode = STATUS_SUCCESS;
		// ヘッダ処理
		if (param->raw == 0) {
			header = malloc(1 * 1024 * 1024);
			if (header == NULL) {
				retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
				break;
			}
			if (param->split == 0) {
				data_size = (param->sampling * (param->bitwidth / 8)) * param->ch;
				data_size *= 10;	// 10秒
			} else {
				data_size = (param->sampling * (param->bitwidth / 8)) * 1;
				data_size *= 10;	// 10秒
			}
			buf1 = malloc(data_size);
			if (buf1 == NULL) {
				retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
				break;
			}
			if (param->ch >= 2) {
				buf2 = malloc(data_size);
				if (buf2 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->ch >= 3) {
				buf3 = malloc(data_size);
				if (buf3 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->ch >= 4) {
				buf4 = malloc(data_size);
				if (buf4 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->ch >= 5) {
				buf5 = malloc(data_size);
				if (buf5 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->ch >= 6) {
				buf6 = malloc(data_size);
				if (buf6 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->rf64 == 0) {
				retCode = PLG_MakeHeaderWAV(inFmt,outFmt,header,1 * 1024 * 1024,&header_size);
				if (retCode != STATUS_SUCCESS) {
					break;
				}
				if (param->ch >= 3 && param->split == 0) {
					chMask	= 0;
					chMask |= SPEAKER_FRONT_LEFT;
					chMask |= SPEAKER_FRONT_RIGHT;
					if (param->chC == 1) {
						chMask |= SPEAKER_FRONT_CENTER;
					}
					if (param->chS == 1) {
						chMask |= SPEAKER_BACK_LEFT;
						chMask |= SPEAKER_BACK_RIGHT;
					}
					if (param->chLFE == 1) {
						chMask |= SPEAKER_LOW_FREQUENCY;
					}
					(*(DWORD *)(&header[40])) = chMask;
				}
			} else {
				retCode = PLG_MakeHeaderRF64(inFmt,outFmt,header,1 * 1024 * 1024,&header_size);
				if (retCode != STATUS_SUCCESS) {
					break;
				}
			}
			if (param->split == 0) {
				ws = fio_write(header,1,header_size,param->fp_w[0]);
				if (param->fp_w[0]->error || ws != header_size) {
					retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
					break;
				}
				fio_get_datasize(param->fp_w[0],&data_start);
			} else {
				for (ccc = 0;ccc < param->ch;ccc++) {
					ws = fio_write(header,1,header_size,param->fp_w[ccc]);
					if (param->fp_w[ccc]->error || ws != header_size) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
					fio_get_datasize(param->fp_w[ccc],&data_start);
				}
				if (retCode != STATUS_SUCCESS) {
					break;
				}
			}
		}
		p_fp1 = fp_r1;
		p_fp2 = fp_r2;
		if (param->ch == 3) {
			if (param->chC == 1) {
				p_fp3 = fp_r3;
			} else if (param->chLFE == 1) {
				p_fp3 = fp_r6;
			}
		}
		if (param->ch == 4) {
			if (param->chS == 1) {
				p_fp3 = fp_r4;
				p_fp4 = fp_r5;
			}
			if (param->chC == 1 && param->chLFE == 1) {
				p_fp3 = fp_r3;
				p_fp4 = fp_r6;
			}
		}
		if (param->ch == 5) {
			if (param->chC == 1 && param->chS == 1) {
				p_fp3 = fp_r3;
				p_fp4 = fp_r4;
				p_fp5 = fp_r5;
			}
			if (param->chLFE == 1 && param->chS == 1) {
				p_fp3 = fp_r6;
				p_fp4 = fp_r4;
				p_fp5 = fp_r5;
			}
		}
		if (param->ch == 6) {
			p_fp3 = fp_r3;
			p_fp4 = fp_r6;
			p_fp5 = fp_r4;
			p_fp6 = fp_r5;
		}
		per = -1;
		for (i = *startSample;i < nSample;i += param->sampling * 10) {
			persent = ((double)i / nSample);
			persent *= 100;
			persent = (int)persent;
			if (persent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;

			// 10秒ずつ処理する
			outSample = param->sampling * 10;
			if (i + outSample > nSample) {
				outSample = nSample - i;
			}
			if (buf1) memset(buf1,0,data_size);
			if (buf2) memset(buf2,0,data_size);
			if (buf3) memset(buf3,0,data_size);
			if (buf4) memset(buf4,0,data_size);
			if (buf5) memset(buf5,0,data_size);
			if (buf6) memset(buf6,0,data_size);

			// データ順
			// Left,Right,Center,LFE,Back Left,Back Right

			#pragma omp parallel
			{
				#pragma omp sections
				{
					#pragma omp section
					{
						retCode1 = Normalize_Mx(param->ch,0,param->bitwidth,p_fp1,outSample,buf1,param);
					}
					#pragma omp section
					{
						if (param->ch >= 2) {
							retCode2 = Normalize_Mx(param->ch,1,param->bitwidth,p_fp2,outSample,buf2,param);
						}
					}
					#pragma omp section
					{
						if (param->ch >= 3) {
							retCode3 = Normalize_Mx(param->ch,2,param->bitwidth,p_fp3,outSample,buf3,param);
						}
					}
					#pragma omp section
					{
						if (param->ch >= 4) {
							retCode4 = Normalize_Mx(param->ch,3,param->bitwidth,p_fp4,outSample,buf4,param);
						}
					}
					#pragma omp section
					{
						if (param->ch >= 5) {
							retCode5 = Normalize_Mx(param->ch,4,param->bitwidth,p_fp5,outSample,buf5,param);
						}
					}
					#pragma omp section
					{
						if (param->ch >= 6) {
							retCode6 = Normalize_Mx(param->ch,5,param->bitwidth,p_fp6,outSample,buf6,param);
						}
					}
				}
			}
			retCode = retCode1 | retCode2 | retCode3 | retCode4 | retCode5 | retCode6;
			if (retCode != STATUS_SUCCESS) {
				break;
			}
			if (param->split == 0) {
				if (param->ch == 2) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j];
					}
				} else if (param->ch == 3) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j] | buf3[j];
					}
				} else if (param->ch == 4) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j] | buf3[j] | buf4[j];
					}
				} else if (param->ch == 5) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j] | buf3[j] | buf4[j] | buf5[j];
					}
				} else if (param->ch == 6) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j] | buf3[j] | buf4[j] | buf5[j] | buf6[j];
					}
				}
				if (outSample == param->sampling * 10) {
					ws = fio_write(buf1,1,data_size,param->fp_w[0]);
					if (param->fp_w[0]->error || ws != data_size) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
				} else {
					write_size = (outSample * (param->bitwidth / 8)) * param->ch;
					ws = fio_write(buf1,1,write_size,param->fp_w[0]);
					if (param->fp_w[0]->error || ws != write_size) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
				}
			} else {
				if (outSample == param->sampling * 10) {
					BYTE *pb[6];
					pb[0] = buf1;
					pb[1] = buf2;
					pb[2] = buf3;
					pb[3] = buf4;
					pb[4] = buf5;
					pb[5] = buf6;
					for (ccc = 0;ccc < param->ch;ccc++) {
						ws = fio_write(pb[ccc],1,data_size,param->fp_w[ccc]);
						if (param->fp_w[ccc]->error || ws != data_size) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					}
					if (retCode != STATUS_SUCCESS) {
						break;
					}
				} else {
					BYTE *pb[6];
					pb[0] = buf1;
					pb[1] = buf2;
					pb[2] = buf3;
					pb[3] = buf4;
					pb[4] = buf5;
					pb[5] = buf6;
					for (ccc = 0;ccc < param->ch;ccc++) {
						write_size = (outSample * (param->bitwidth / 8)) * 1;
						ws = fio_write(pb[ccc],1,write_size,param->fp_w[ccc]);
						if (param->fp_w[ccc]->error || ws != write_size) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					}
					if (retCode != STATUS_SUCCESS) {
						break;
					}
				}
			}
			if (param->rf64 == 0) {
				fio_get_datasize(param->fp_w[0],&file_size);
				if (file_size > (fio_size)1 * 1024 * 1024 * 1024 && outSample == param->sampling * 10 && (nSample - i) >= param->sampling * 20) {
					// データサイズが大きいので分割する(rf64以外)
					*startSample = i + outSample;
					div_flag = 1;
					break;
				}
			}
		}

		if (retCode != STATUS_SUCCESS) {
			break;
		}
		fio_get_datasize(param->fp_w[0],&data_end);
		bwf_enable = 0;
		if (param->bwf) {
			bwf_enable = 1;
		}
		if (param->raw == 0) {
			chunk_data = NULL;
			for (nChunk = 1;;) {
				PLG_GetExtraChunk(param->fromfile,nChunk,&chunk_data,&chunk_size);
if (1) {
char s[256];
sprintf(s,"%s:chunk_size:%d",param->fromfile,chunk_size);
PRINT_LOG(s);
}
				if (chunk_size == 0) {
					break;
				}
				if (chunk_size > 0 && chunk_data[0] == 'b' && chunk_data[1] == 'e' && chunk_data[2] == 'x' && chunk_data[3] == 't') {
					bwf_enable = 0;
					if (param->bwf) {
						int alloc_flag = 0;
						bwf_size = chunk_data[7];bwf_size <<= 8;
						bwf_size |= chunk_data[6];bwf_size <<= 8;
						bwf_size |= chunk_data[5];bwf_size <<= 8;
						bwf_size |= chunk_data[4];
						if (bwf_size <= sizeof (BROADCAST_EXT) + 128) {
							alloc_flag = 1;
						} else {
							bext = (BROADCAST_EXT *)&chunk_data[8];
							if (chunk_size < sizeof (BROADCAST_EXT) + strlen(bext->codingHistory) + 128) {
								alloc_flag = 1;
							}
						}
						if (alloc_flag == 1) {
							bwf_size += 256;
							bext = (BROADCAST_EXT *)malloc(bwf_size);
						}
						if (bext) {
							if (alloc_flag == 1) {
								memset(bext,0,bwf_size);
								memcpy(bext,chunk_data + 8,chunk_size);
							}
							UpdateBext(bext,inFmt,outFmt,param,bwf_size);
							if (param->split == 0) {
								ws = fio_write("bext",1,4,param->fp_w[0]);
								if (param->fp_w[0]->error || ws != 4) {
									retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
									break;
								}
								ws = fio_write(&bwf_size,1,4,param->fp_w[0]);
								if (param->fp_w[0]->error || ws != 4) {
									retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
									break;
								}
								ws = fio_write(bext,1,bwf_size,param->fp_w[0]);
								if (param->fp_w[0]->error || ws != bwf_size) {
									retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
									break;
								}
							} else {
								for (ccc = 0;ccc < param->ch;ccc++) {
									ws = fio_write("bext",1,4,param->fp_w[ccc]);
									if (param->fp_w[ccc]->error || ws != 4) {
										retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
										break;
									}
									ws = fio_write(&bwf_size,1,4,param->fp_w[ccc]);
									if (param->fp_w[ccc]->error || ws != 4) {
										retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
										break;
									}
									ws = fio_write(bext,1,bwf_size,param->fp_w[ccc]);
									if (param->fp_w[ccc]->error || ws != bwf_size) {
										retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
										break;
									}
								}
								if (retCode != STATUS_SUCCESS) {
									break;
								}
							}
							if (alloc_flag) {
								free(bext);
							}
						}
					}
				} else {
					if (chunk_size > 0) {
						if (param->split == 0) {
							int sz;
							sz = chunk_size;
							if (sz & 1) sz++;
							ws = fio_write(chunk_data,1,sz,param->fp_w[0]);
						} else {
							int sz;
							sz = chunk_size;
							if (sz & 1) sz++;
							for (ccc = 0;ccc < param->ch;ccc++) {
								ws = fio_write(chunk_data,1,sz,param->fp_w[ccc]);
								if (param->fp_w[ccc]->error || ws != sz) {
									retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
									break;
								}
							}
							if (retCode != STATUS_SUCCESS) {
								break;
							}
						}
					}
				}
				nChunk++;
			}
			if (bwf_enable) {
				bwf_chunk = malloc(512 + 256);
				if (bwf_chunk) {
					memset(bwf_chunk,0,512 + 256);
					bwf_chunk[0] = 'b';
					bwf_chunk[1] = 'e';
					bwf_chunk[2] = 'x';
					bwf_chunk[3] = 't';
					bwf_chunk[4] = (BYTE)((512 + 256) - 8);
					bwf_chunk[5] = (BYTE)(((512 + 256) - 8) >> 8);
					bwf_chunk[6] = (BYTE)(((512 + 256) - 8) >> 16);
					bwf_chunk[7] = (BYTE)(((512 + 256) - 8) >> 24);
					bext = (BROADCAST_EXT *)&bwf_chunk[8];
					UpdateBext(bext,inFmt,outFmt,param,512 + 256);
					if (param->split == 0) {
						ws = fio_write(bwf_chunk,1,512 + 256,param->fp_w[0]);
						if (param->fp_w[0]->error || ws != 512 + 256) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					} else {
						for (ccc = 0;ccc < param->ch;ccc++) {
							ws = fio_write(bwf_chunk,1,512 + 256,param->fp_w[ccc]);
							if (param->fp_w[ccc]->error || ws != 512 + 256) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						}
						if (retCode != STATUS_SUCCESS) {
							break;
						}
					}
					free(bwf_chunk);
					bwf_chunk = NULL;
				}
			}
			if (param->bwf && (div_flag == 1 || *pCount > 1)) {
				// bwf の link チャンク
				char *link,*wk_str;
				long link_size;
				link = (char *)malloc(strlen(param->tofile) + 128);
				wk_str = malloc(_MAX_PATH + 128);
				if (link != NULL && wk_str != NULL) {
					link[0] = '\0';
					strcat(link,link_start);
					if (*pCount == 1) {
						sprintf(wk_str,link_file,"actual",(int)*pCount,param->tofile);
					} else {
						sprintf(wk_str,link_file,"other",(int)*pCount,param->tofile);
					}
					strcat(link,wk_str);
					strcat(link,link_end);
					if (param->split == 0) {
						ws = fio_write("link",1,4,param->fp_w[0]);
						if (param->fp_w[0]->error || ws != 4) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					} else {
						for (ccc = 0;ccc < param->ch;ccc++) {
							ws = fio_write("link",1,4,param->fp_w[ccc]);
							if (param->fp_w[ccc]->error || ws != 4) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						}
						if (retCode != STATUS_SUCCESS) {
							break;
						}
					}
					link_size = strlen(link);
					if (param->split == 0) {
						ws = fio_write(&link_size,1,4,param->fp_w[0]);
						if (param->fp_w[0]->error || ws != 4) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					} else {
						for (ccc = 0;ccc < param->ch;ccc++) {
							ws = fio_write(&link_size,1,4,param->fp_w[ccc]);
							if (param->fp_w[ccc]->error || ws != 4) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						}
					}
					if (param->split == 0) {
						ws = fio_write(link,1,link_size,param->fp_w[0]);
						if (param->fp_w[0]->error || ws != link_size) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					} else {
						for (ccc = 0;ccc < param->ch;ccc++) {
							ws = fio_write(link,1,link_size,param->fp_w[ccc]);
							if (param->fp_w[ccc]->error || ws != link_size) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						}
						if (retCode != STATUS_SUCCESS) {
							break;
						}
					}
					free(link);
					free(wk_str);
					if (link_size & 0x01) {
						// padding
						if (param->split == 0) {
							ws = fio_write("",1,1,param->fp_w[0]);
							if (param->fp_w[0]->error || ws != 1) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						} else {
							for (ccc = 0;ccc < param->ch;ccc++) {
								ws = fio_write("",1,1,param->fp_w[ccc]);
								if (param->fp_w[ccc]->error || ws != 1) {
									retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
									break;
								}
							}
						}
					}
				}
			}
			if (param->split == 0) {
				fio_get_datasize(param->fp_w[0],&data_size);
				fio_flush(param->fp_w[0]);
				if (param->fp_w[0]->error) {
					retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
					break;
				}
				fio_rewind(param->fp_w[0]);
				if (param->rf64 == 0) {
					retCode = PLG_UpdateHeaderWAV(outFmt,data_size,data_end - data_start,header,header_size);
				} else {
					retCode = PLG_UpdateHeaderRF64(outFmt,data_size,data_end - data_start,header,header_size);
				}
if (1) {
char s[256];
sprintf(s,"outSize:%lld\n",data_end - data_start);
PRINT_LOG(s);
}
				ws = fio_write(header,1,header_size,param->fp_w[0]);
				if (param->fp_w[0]->error || ws != header_size) {
					retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
					break;
				}
			} else {
				fio_get_datasize(param->fp_w[0],&data_size);
				for (ccc = 0;ccc < param->ch;ccc++) {
					fio_flush(param->fp_w[ccc]);
					if (param->fp_w[ccc]->error) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
					fio_rewind(param->fp_w[ccc]);
					if (param->rf64 == 0) {
						retCode = PLG_UpdateHeaderWAV(outFmt,data_size,data_end - data_start,header,header_size);
					} else {
						retCode = PLG_UpdateHeaderRF64(outFmt,data_size,data_end - data_start,header,header_size);
					}
					ws = fio_write(header,1,header_size,param->fp_w[ccc]);
					if (param->fp_w[ccc]->error || ws != header_size) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
				}
			}
		}
	} while (0);
	if (div_flag == 0) {
		*pCount = 0;
	}
	fprintf(stdout,"%d%%\n",100);
	fflush(stdout);

	return retCode;
}
//---------------------------------------------------------------------------
// Function   : Normalize_Mx
// Description: ノーマライズ処理
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_Mx(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	int retCode;

	// 音声データ出力
	switch (param->mode) {
		case 0:
			retCode = Normalize_M0(nCh,ch,bit,fp_r,nSample,buffer,param);
			break;
		case 1:
			retCode = Normalize_M1(nCh,ch,bit,fp_r,nSample,buffer,param);
			break;
		case 2:
			retCode = Normalize_M2(nCh,ch,bit,fp_r,nSample,buffer,param);
			break;
		case 3:
			retCode = Normalize_M3(nCh,ch,bit,fp_r,nSample,buffer,param);
			break;
		default:
			retCode = -1;
			break;
	}
	return retCode;
}
//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理(カット)
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_M0(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	SSIZE s;
	DWORD i,w_off;
	int next;
	double nx;
	short *s_os;
	float *f_os;
	double *d_os;
	fio_size rs;
	char p[50];
	if (bit == 16) {
		next = 2 * nCh;
		w_off = 2 * ch;
		if (param->split) {
			next = 2;
			w_off = 0;
		}
	} else if (bit == 24) {
		next = 3 * nCh;
		w_off = 3 * ch;
		if (param->split) {
			next = 3;
			w_off = 0;
		}
	} else if (bit == 32) {
		next = 4 * nCh;
		w_off = 4 * ch;
		if (param->split) {
			next = 4;
			w_off = 0;
		}
	} else if (bit == 64) {
		next = 8 * nCh;
		w_off = 8 * ch;
		if (param->split) {
			next = 8;
			w_off = 0;
		}
	}

	if (ch == 0) {
		nx = param->nx1;
	} else if (ch == 1) {
		nx = param->nx2;
	} else if (ch == 2) {
		nx = param->nx3;
	} else if (ch == 3) {
		nx = param->nx4;
	} else if (ch == 4) {
		nx = param->nx5;
	} else if (ch == 5) {
		nx = param->nx6;
	}
	for (i = 0;i < nSample;i++) {
		rs = fio_read(&s,sizeof (SSIZE),1,fp_r);
		if (fp_r->error || rs != 1) {
			param->errLine = __LINE__;
			return STATUS_FILE_READ_ERR;
		}
		s = CLIP_NX(s,nx);
		
		if (bit == 16) {
			s = CLIP_MAX(s);
			s >>= (64 - 16);
			s_os = (short *)&buffer[w_off];
			*s_os = (short)s;
			w_off += next;
		} else if (bit == 24) {
			s = CLIP_MAX(s);
			s >>= (64 - 24);
			buffer[w_off + 0] = (unsigned char)s;
			buffer[w_off + 1] = (unsigned char)(s >> 8);
			buffer[w_off + 2] = (unsigned char)(s >> 16);
			w_off += next;
		} else if (bit == 32) {
			f_os = (float *)&buffer[w_off];
			s = CLIP_MAX(s);
			*f_os = (float)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		} else if (bit == 64) {
			d_os = (double *)&buffer[w_off];
			s = CLIP_MAX(s);
			*d_os = (double)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		}
	}
	return STATUS_SUCCESS;
}
//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理(ディザ)
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_M1(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	SSIZE s;
	DWORD i;
	int ignore_s;
	double nx;
	double noise;
	short *s_os;
	float *f_os;
	double *d_os;
	fio_size rs;
	int next;
	int w_off;
	
	if (bit == 16) {
		next = 2 * nCh;
		w_off = 2 * ch;
		if (param->split) {
			next = 2;
			w_off = 0;
		}
	} else if (bit == 24) {
		next = 3 * nCh;
		w_off = 3 * ch;
		if (param->split) {
			next = 3;
			w_off = 0;
		}
	} else if (bit == 32) {
		next = 4 * nCh;
		w_off = 4 * ch;
		if (param->split) {
			next = 4;
			w_off = 0;
		}
	} else if (bit == 64) {
		next = 8 * nCh;
		w_off = 8 * ch;
		if (param->split) {
			next = 8;
			w_off = 0;
		}
	}

	if (ch == 0) {
		nx = param->nx1;
	} else if (ch == 1) {
		nx = param->nx2;
	} else if (ch == 2) {
		nx = param->nx3;
	} else if (ch == 3) {
		nx = param->nx4;
	} else if (ch == 4) {
		nx = param->nx5;
	} else if (ch == 5) {
		nx = param->nx6;
	}

	for (i = 0;i < nSample;i++) {
		rs = fio_read(&s,sizeof (SSIZE),1,fp_r);
		if (fp_r->error || rs != 1) {
			param->errLine = __LINE__;
			return STATUS_FILE_READ_ERR;
		}

		s = CLIP_NX(s,nx);

		noise = normalNoise();
		if (bit == 16) {
			if (param->ditherLv > 0) {
				noise *= (0x10000000000 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-16));
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX(s);

			s >>= (64 - 16);
			s_os = (short *)&buffer[w_off];
			*s_os = (short)s;
			w_off += next;

		} else if (bit == 24) {
			if (param->ditherLv > 0) {
				noise *= (0x100000000 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-24));
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX(s);

			s >>= (64 - 24);
			buffer[w_off + 0] = (unsigned char)s;
			buffer[w_off + 1] = (unsigned char)(s >> 8);
			buffer[w_off + 2] = (unsigned char)(s >> 16);
			w_off += next;
		} else if (bit == 32) {
			if (param->ditherLv > 0) {
				noise *= (0x1000000 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-32));
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX(s);

			f_os = (float *)&buffer[w_off];
			*f_os = (float)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		} else if (bit == 64) {
			if (param->ditherLv > 0) {
				noise *= (0x100 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-48));
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX(s);

			d_os = (double *)&buffer[w_off];
			*d_os = (double)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		}
	}

	return STATUS_SUCCESS;
}
//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理(誤差累積)
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_M2(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	SSIZE s;
	SSIZE old_s;
	SSIZE sum,sum_2nd,sum_3rd,val;
	DWORD i;
	int ignore_s;
	double nx;
	double noise;
	short *s_os;
	float *f_os;
	double *d_os;
	fio_size rs;
	int next;
	int w_off;
	char p[100];


	if (bit == 16) {
		next = 2 * nCh;
		w_off = 2 * ch;
		if (param->split) {
			next = 2;
			w_off = 0;
		}
	} else if (bit == 24) {
		next = 3 * nCh;
		w_off = 3 * ch;
		if (param->split) {
			next = 3;
			w_off = 0;
		}
	} else if (bit == 32) {
		next = 4 * nCh;
		w_off = 4 * ch;
		if (param->split) {
			next = 4;
			w_off = 0;
		}
	} else if (bit == 64) {
		next = 8 * nCh;
		w_off = 8 * ch;
		if (param->split) {
			next = 8;
			w_off = 0;
		}
	}

	if (ch == 0) {
		nx = param->nx1;
		sum = param->sum1;
		sum_2nd = param->sum1_2nd;
		sum_3rd = param->sum1_3rd;
		old_s	= param->old_s1;
	} else if (ch == 1) {
		nx = param->nx2;
		sum = param->sum2;
		sum_2nd = param->sum2_2nd;
		sum_3rd = param->sum2_3rd;
		old_s	= param->old_s2;
	} else if (ch == 2) {
		nx = param->nx3;
		sum = param->sum3;
		sum_2nd = param->sum3_2nd;
		sum_3rd = param->sum3_3rd;
		old_s	= param->old_s3;
	} else if (ch == 3) {
		nx = param->nx4;
		sum = param->sum4;
		sum_2nd = param->sum4_2nd;
		sum_3rd = param->sum4_3rd;
		old_s	= param->old_s4;
	} else if (ch == 4) {
		nx = param->nx5;
		sum = param->sum5;
		sum_2nd = param->sum5_2nd;
		sum_3rd = param->sum5_3rd;
		old_s	= param->old_s5;
	} else if (ch == 5) {
		nx = param->nx6;
		sum = param->sum6;
		sum_2nd = param->sum6_2nd;
		sum_3rd = param->sum6_3rd;
		old_s	= param->old_s6;
	}

	for (i = 0;i < nSample;i++) {
		rs = fio_read(&s,sizeof (SSIZE),1,fp_r);
		if (fp_r->error || rs != 1) {
			param->errLine = __LINE__;
			return STATUS_FILE_READ_ERR;
		}

		s = CLIP_NX(s,nx);

		noise = normalNoise();
		if (bit == 16) {
			// 四捨五入

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-16));
			
			s -= old_s;
			sum += s;

			s = sum >> ((64-8) - 16);
			s <<= ((64-8) - 16);
			s -= old_s;
			sum_2nd += s;
			s = sum_2nd >> ((64-8) - 16);

			s <<= ((64-8) - 16);
			s -= old_s;
			sum_3rd += s;
			s = sum_3rd >> ((64-8) - 16);

			old_s = s << ((64-8) - 16);

			if (param->ditherLv > 0) {
				noise *= (1 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX_N(s,16);
			s_os = (short *)&buffer[w_off];
			*s_os = (short)s;
			w_off += next;
		} else if (bit == 24) {

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-24));
			s -= old_s;
			sum += s;

			s = sum >> ((64-8) - 24);
			s <<= ((64-8) - 24);
			s -= old_s;
			sum_2nd += s;
			s = sum_2nd >> ((64-8) - 24);

			s <<= ((64-8) - 24);
			s -= old_s;
			sum_3rd += s;
			s = sum_3rd >> ((64-8) - 24);

			old_s = s << ((64-8) - 24);

			if (param->ditherLv > 0) {
				noise *= (1 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX_N(s,24);
			buffer[w_off + 0] = (unsigned char)s;
			buffer[w_off + 1] = (unsigned char)(s >> 8);
			buffer[w_off + 2] = (unsigned char)(s >> 16);
			w_off += next;
		} else if (bit == 32) {
			s = CLIP_ADD(s,ROUND_NBIT((64-8)-32));

			// 誤差の累積
			sum += s % 0x1000000;

			// 誤差の判定
			val = sum / 0x1000000;
			if (ignore_s == 0 && val > 0) {
				s += 0x1000000;
				sum %= 0x1000000;
			} else if (ignore_s == 0 && val < 0) {
				s -= 0x1000000;
				sum %= 0x1000000;
			}

			// 2次 誤差の累積
			sum_2nd += s % 0x1000000;
			// 誤差の判定
			val = sum_2nd / 0x1000000;
			if (val > 0) {
				s += 0x1000000;
				sum_2nd %= 0x1000000;
			} else if (val < 0) {
				s -= 0x1000000;
				sum_2nd %= 0x1000000;
			}

			// 32bit 化する
			s >>= ((64-8) - 32);
			if (param->ditherLv > 0) {
				noise *= (1 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX_N(s,32);
			s <<= (64 - 32);
			f_os = (float *)&buffer[w_off];
			*f_os = (float)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		} else if (bit == 64) {
			// 48bit
			s = CLIP_ADD(s,ROUND_NBIT((64-8)-48));

			// 誤差の累積
			sum += s % 0x10000;

			// 誤差の判定
			val = sum / 0x10000;
			if (val > 0) {
				s += 0x10000;
				sum %= 0x10000;
			} else if (val < 0) {
				s -= 0x10000;
				sum %= 0x10000;
			}

			// 2次 誤差の累積
			sum_2nd += s % 0x10000;
			// 誤差の判定
			val = sum_2nd / 0x10000;
			if (val > 0) {
				s += 0x10000;
				sum_2nd %= 0x10000;
			} else if (val < 0) {
				s -= 0x10000;
				sum_2nd %= 0x10000;
			}

			// 64bit 化する
			s >>= ((64-8) - 48);
			if (param->ditherLv > 0) {
				noise *= (1 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX_N(s,48);
			s <<= (64 - 48);
			d_os = (double *)&buffer[w_off];
			*d_os = (double)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		}
	}
	if (ch == 0) {
		param->sum1 = sum;
		param->sum1_2nd = sum_2nd;
		param->sum1_3rd = sum_3rd;
		param->old_s1 = old_s;
	} else if (ch == 1) {
		param->sum2 = sum;
		param->sum2_2nd = sum_2nd;
		param->sum2_3rd = sum_3rd;
		param->old_s2 = old_s;
	} else if (ch == 2) {
		param->sum3 = sum;
		param->sum3_2nd = sum_2nd;
		param->sum3_3rd = sum_3rd;
		param->old_s3 = old_s;
	} else if (ch == 3) {
		param->sum4 = sum;
		param->sum4_2nd = sum_2nd;
		param->sum4_3rd = sum_3rd;
		param->old_s4 = old_s;
	} else if (ch == 4) {
		param->sum5 = sum;
		param->sum5_2nd = sum_2nd;
		param->sum5_3rd = sum_3rd;
		param->old_s5 = old_s;
	} else if (ch == 5) {
		param->sum6 = sum;
		param->sum6_2nd = sum_2nd;
		param->sum6_3rd = sum_3rd;
		param->old_s6 = old_s;
	}

	return STATUS_SUCCESS;
}
//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理(誤差拡散)
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_M3(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	SSIZE s;
	SSIZE ss[3];
	SSIZE sd[3];
	SSIZE val;
	SSIZE *buf;
	DWORD pastSample;
	long i,remain;
	int ignore_s;
	double nx;
	double noise;
	short *s_os;
	float *f_os;
	double *d_os;
	fio_size rs;
	int next;
	int w_off;

	if (bit == 16) {
		next = 2 * nCh;
		w_off = 2 * ch;
		if (param->split) {
			next = 2;
			w_off = 0;
		}
	} else if (bit == 24) {
		next = 3 * nCh;
		w_off = 3 * ch;
		if (param->split) {
			next = 3;
			w_off = 0;
		}
	} else if (bit == 32) {
		next = 4 * nCh;
		w_off = 4 * ch;
		if (param->split) {
			next = 4;
			w_off = 0;
		}
	} else if (bit == 64) {
		next = 8 * nCh;
		w_off = 8 * ch;
		if (param->split) {
			next = 8;
			w_off = 0;
		}
	}

	if (ch == 0) {
		nx = param->nx1;
	} else if (ch == 1) {
		nx = param->nx2;
	} else if (ch == 2) {
		nx = param->nx3;
	} else if (ch == 3) {
		nx = param->nx4;
	} else if (ch == 4) {
		nx = param->nx5;
	} else if (ch == 5) {
		nx = param->nx6;
	}

	ss[0] = ss[1] = ss[2] = 0;
	sd[0] = sd[1] = sd[2] = 0;

	buf = malloc(nSample * sizeof (SSIZE));
	if (buf == NULL) {
		param->errLine = __LINE__;
		return STATUS_MEM_ALLOC_ERR;
	}
	
	pastSample = 0;
	do {
		remain = fio_read(buf,sizeof (SSIZE),nSample,fp_r);
		for (i = 1;i < remain + 1 && pastSample < nSample;i++,pastSample++) {
			noise = normalNoise();

			ss[0] = buf[i-1];
			if (i < remain) {
				ss[1] = buf[i];
			} else {
				ss[1] = ss[0];
			}
			if (i + 1 < remain) {
				ss[2] = buf[i+1];
			} else {
				ss[2] = ss[0];
			}

			ss[0] = CLIP_NX(ss[0],nx);
			ss[1] = CLIP_NX(ss[1],nx);
			ss[2] = CLIP_NX(ss[2],nx);

			if (bit == 16) {
				ss[1] += sd[1];
				s = ss[1];
				val = ss[1] % 0x10000000000;
				ss[1] = ss[1] - val;
				if (val > 0) {
					if (val > 0x7FFFFFFFFF) {
						ss[1] += 0x10000000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFFFFFF) {
						ss[1] -= 0x10000000000;
					}
				}
				sd[1]  = ss[1] - s;
				sd[0] += sd[1] / 2;
				sd[2] += sd[1] / 2;
				sd[1] = 0;
				ss[0] += sd[0];
				s = ss[0];
				val = ss[0] % 0x10000000000;
				ss[0] = ss[0] - val;
				if (val > 0) {
					if (val > 0x7FFFFFFFFF) {
						ss[0] += 0x10000000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFFFFFF) {
						ss[0] -= 0x10000000000;
					}
				}
				sd[0] = ss[0] - s;
				sd[1] = sd[0] / 2;

				sd[0] = sd[1];
				sd[1] = sd[2];
				sd[2] = 0;
				s = ss[0];
				s >>= ((64-8) - 16);

				if (param->ditherLv > 0) {
					noise *= (1 << (param->ditherLv - 1));
				} else {
					noise = 0;
				}
				s = CLIP_ADD(s,(SSIZE)noise);
				s = CLIP_MAX_N(s,16);
				s_os = (short *)&buffer[w_off];
				*s_os = (short)s;
				w_off += next;
			} else if (bit == 24) {
				ss[1] += sd[1];
				s = ss[1];
				val = ss[1] % 0x100000000;
				ss[1] = ss[1] - val;
				if (val > 0) {
					if (val > 0x7FFFFFFF) {
						ss[1] += 0x100000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFFFF) {
						ss[1] -= 0x100000000;
					}
				}
				sd[1]  = ss[1] - s;
				sd[0] += sd[1] / 2;
				sd[2] += sd[1] / 2;
				sd[1] = 0;
				ss[0] += sd[0];
				s = ss[0];
				val = ss[0] % 0x100000000;
				ss[0] = ss[0] - val;
				if (val > 0) {
					if (val > 0x7FFFFFFF) {
						ss[0] += 0x100000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFFFF) {
						ss[0] -= 0x100000000;
					}
				}
				sd[0] = ss[0] - s;
				sd[1] = sd[0] / 2;

				sd[0] = sd[1];
				sd[1] = sd[2];
				sd[2] = 0;
				s = ss[0];
				s >>= ((64-8) - 24);

				if (param->ditherLv > 0) {
					noise *= (1 << (param->ditherLv - 1));
				} else {
					noise = 0;
				}
				s = CLIP_ADD(s,(SSIZE)noise);
				s = CLIP_MAX_N(s,24);
				buffer[w_off + 0] = (unsigned char)s;
				buffer[w_off + 1] = (unsigned char)(s >> 8);
				buffer[w_off + 2] = (unsigned char)(s >> 16);
				w_off += next;
			} else if (bit == 32) {
				ss[1] += sd[1];
				s = ss[1];
				val = ss[1] % 0x1000000;
				ss[1] = ss[1] - val;
				if (val > 0) {
					if (val > 0x7FFFFF) {
						ss[1] += 0x1000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFF) {
						ss[1] -= 0x1000000;
					}
				}
				sd[1]  = ss[1] - s;
				sd[0] += sd[1] / 2;
				sd[2] += sd[1] / 2;
				sd[1] = 0;
				ss[0] += sd[0];
				s = ss[0];
				val = ss[0] % 0x1000000;
				ss[0] = ss[0] - val;
				if (val > 0) {
					if (val > 0x7FFFFF) {
						ss[0] += 0x1000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFF) {
						ss[0] -= 0x1000000;
					}
				}
				sd[0] = ss[0] - s;
				sd[1] = sd[0] / 2;

				sd[0] = sd[1];
				sd[1] = sd[2];
				sd[2] = 0;
				s = ss[0];
				s >>= ((64-8) - 32);

				if (param->ditherLv > 0) {
					noise *= (1 << (param->ditherLv - 1));
				} else {
					noise = 0;
				}
				s = CLIP_ADD(s,(SSIZE)noise);
				s = CLIP_MAX_N(s,32);
				s <<= (64 - 32);
				f_os = (float *)&buffer[w_off];
				*f_os = (float)s / 0x7FFFFFFFFFFFFFFF;
				w_off += next;
			} else if (bit == 64) {
				ss[1] += sd[1];
				s = ss[1];
				val = ss[1] % 0x100;
				ss[1] = ss[1] - val;
				if (val > 0) {
					if (val > 0x7F) {
						ss[1] += 0x100;
					}
				} else {
					if ((val * -1) > 0x7F) {
						ss[1] -= 0x100;
					}
				}
				sd[1]  = ss[1] - s;
				sd[0] += sd[1] / 2;
				sd[2] += sd[1] / 2;
				sd[1] = 0;
				ss[0] += sd[0];
				s = ss[0];
				val = ss[0] % 0x100;
				ss[0] = ss[0] - val;
				if (val > 0) {
					if (val > 0x7F) {
						ss[0] += 0x100;
					}
				} else {
					if ((val * -1) > 0x7F) {
						ss[0] -= 0x100;
					}
				}
				sd[0] = ss[0] - s;
				sd[1] = sd[0] / 2;

				sd[0] = sd[1];
				sd[1] = sd[2];
				sd[2] = 0;
				s = ss[0];
				s >>= ((64-8) - 48);

				if (param->ditherLv > 0) {
					noise *= (1 << (param->ditherLv - 1));
				} else {
					noise = 0;
				}
				s = CLIP_ADD(s,(SSIZE)noise);
				s = CLIP_MAX_N(s,48);
				s <<= (64 - 48);
				d_os = (double *)&buffer[w_off];
				*d_os = (double)s / 0x7FFFFFFFFFFFFFFF;
				w_off += next;
			}
		}
	} while (pastSample < nSample);
	free(buf);

	return STATUS_SUCCESS;
}

//---------------------------------------------------------------------------
// Function   : Normalize_DSD
// Description: DSD用ノーマライズ処理
// ---
//	inFmt	:入力ファイル音声形式情報
//	outFmt	:出力ファイル音声形式情報
//	fp_r1	:音声データのFIO構造体
//	fp_r2	:音声データのFIO構造体
//	nSample :処理をするサンプル数のアドレス
//	param	:パラメーター構造体
//
int Normalize_DSD(SOUNDFMT *inFmt,SOUNDFMT *outFmt,FIO *fp_r1,FIO *fp_r2,DWORD nSample,PARAM_INFO *param)
{
	DWORD i,j,outSample;
	DWORD chMask;
	DSF *dsf;
	DSF_FMT *dsf_fmt;
	DSF_DATA *dsf_data;
	int retCode;
	SSIZE *buf1,*buf2;
	BYTE *buf1_4k,*buf2_4k;
	FIO *fp1_w,*fp2_w;
	fio_size write_size,ws,data_size,file_size,rs;
	fio_size data_start,data_end;
	double persent,per;
	SSIZE s;
	double nx1,nx2;
	char *header;
	long hs;

	// 1bit
	SSIZE min,max;
	SSIZE s1,s2;
	SSIZE ss1,ss2;
	SSIZE sss1,sss2;
	SSIZE q1,q2;
	SSIZE d1,d2;

	// 8bit,5bit,3bit
	SSIZE sm1[3],sm2[3];
	SSIZE ssm1[3],ssm2[3];
	SSIZE qm1[3],qm2[3];
	SSIZE dm1[3],dm2[3];
	
	min = pow(2,58) * -1;
	max = pow(2,58);
	s1 = s2 = 0;
	ss1 = ss2 = 0;
	sss1 = sss2 = 0;
	q1 = q2 = 0;

	sm1[0] = sm1[1] = sm1[2] = 0;
	sm2[0] = sm2[1] = sm2[3] = 0;
	ssm1[0] = ssm1[1] = ssm1[2] = 0;
	ssm2[0] = ssm2[1] = ssm2[2] = 0;
	qm1[0] = qm1[1] = qm1[2] = 0;
	qm2[0] = qm2[1] = qm2[2] = 0;
	dm1[0] = dm1[1] = dm1[2] = 0;
	dm2[0] = dm2[1] = dm2[2] = 0;

	nx1 = param->nx1;
	nx2 = param->nx2;

	do {
		retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
		header = malloc(1024);
		buf1 = malloc((4096 * 8) * sizeof (SSIZE));
		buf2 = malloc((4096 * 8) * sizeof (SSIZE));
		buf1_4k = malloc(4096);
		buf2_4k = malloc(4096);
		if (header == NULL || buf1 == NULL || buf2 == NULL || buf1_4k == NULL || buf2_4k == NULL) {
			break;
		}
		memset(header,0,1024);
		
		retCode = PLG_MakeHeaderDSD(inFmt,outFmt,header,1024,&hs);
		if (retCode != STATUS_SUCCESS) {
			param->errLine = __LINE__;
			break;
		}
		ws = fio_write(header,1,hs,param->fp_w[0]);
		if (param->fp_w[0]->error || ws != hs) {
			retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
			break;
		}
		fio_get_datasize(param->fp_w[0],&data_start);

		retCode = STATUS_SUCCESS;
		per = -1;
		data_size = (4096 * 8) * sizeof (SSIZE);
		for (i = 0;i < nSample;i += 4096 * 8) {
			persent = ((double)i / nSample);
			persent *= 100;
			if (persent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;
PRINT_LOG("");
			// 1秒ずつ処理する
			outSample = 4096 * 8;
			if (i + outSample > nSample) {
				outSample = nSample - i;
			}
			if (buf1) memset(buf1,0,data_size);
			if (buf2) memset(buf2,0,data_size);
PRINT_LOG("");
			if (param->ch == 1) {
				rs = fio_read(buf1,sizeof (SSIZE),outSample,fp_r1);
				if (fp_r1->error || rs <= 0) {
					param->errLine = __LINE__;
					return STATUS_FILE_READ_ERR;
				}
				// delta sigma modulation (1bit)
				for (j = 0;j < outSample;j++) {
					s = CLIP_NX(buf1[j],nx1);
					s1  += s - q1;
					ss1 += s1 - q1 * 2;
					
					q1 = (ss1 >= 0) ? max : min;
					buf1[j] = q1;
				}
				nbit2onebit(buf1,buf1_4k,4096 * 8);
				ws = fio_write(buf1_4k,1,4096,param->fp_w[0]);
				if (param->fp_w[0]->error || ws != 4096) {
					retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
					break;
				}
			} else if (param->ch == 2) {
				rs = fio_read(buf1,sizeof (SSIZE),outSample,fp_r1);
				if (fp_r1->error || rs <= 0) {
					param->errLine = __LINE__;
					return STATUS_FILE_READ_ERR;
				}
				// delta sigma modulation (1bit)
				for (j = 0;j < outSample;j++) {
					s = CLIP_NX(buf1[j],nx1);
					s1  += s - q1;
					ss1 += s1 - q1 * 2;
					
					q1 = (ss1 >= 0) ? max : min;
					buf1[j] = q1;
				}
				// 2値化
				nbit2onebit(buf1,buf1_4k,4096 * 8);
				ws = fio_write(buf1_4k,1,4096,param->fp_w[0]);
				if (param->fp_w[0]->error || ws != 4096) {
					retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
					break;
				}

				rs = fio_read(buf2,sizeof (SSIZE),outSample,fp_r2);
				if (fp_r2->error || rs <= 0) {
					param->errLine = __LINE__;
					return STATUS_FILE_READ_ERR;
				}
				// delta sigma modulation (1bit)
				for (j = 0;j < outSample;j++) {
					s = CLIP_NX(buf2[j],nx2);
					s2  += s - q2;
					ss2 += s2 - q2 * 2;
					
					q2 = (ss2 >= 0) ? max : min;
					buf2[j] = q2;
				}
				nbit2onebit(buf2,buf2_4k,4096 * 8);
				ws = fio_write(buf2_4k,1,4096,param->fp_w[0]);
				if (param->fp_w[0]->error || ws != 4096) {
					retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
					break;
				}
			}
		}
		if (retCode != STATUS_SUCCESS) break;

PRINT_LOG("");
		fio_flush(param->fp_w[0]);
		if (param->fp_w[0]->error) {
			retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
			break;
		}
		fio_get_datasize(param->fp_w[0],&data_end);
		fio_rewind(param->fp_w[0]);

		retCode = PLG_UpdateHeaderDSD(outFmt,data_end,data_end - data_start,nSample,header,hs);
		ws = fio_write(header,1,hs,param->fp_w[0]);
		if (param->fp_w[0]->error || ws != hs) {
			retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
			break;
		}
		retCode = STATUS_SUCCESS;
	} while (0);

PRINT_LOG("");
	if (header) free(header);
	if (buf1) free(buf1);
	if (buf2) free(buf2);
	if (buf1_4k) free(buf1_4k);
	if (buf2_4k) free(buf2_4k);

	return retCode;
}
//---------------------------------------------------------------------------
// Function   : nbit2onebit
// Description: 1bit化
//
void nbit2onebit(SSIZE *i_buf,BYTE *o_buf,int size)
{
	int i;
	int bit,ptr_byte;
	if (size > 0) {
		ptr_byte = 0;
		bit  = 0x01;
		for (i = 0;i < size;i++) {
			if ((i_buf[i]) > 0) {
				o_buf[ptr_byte] |= bit;
			} else if ((i_buf[i]) < 0) {
				o_buf[ptr_byte] &= (bit ^ 0xFF);
			}
			bit <<= 1;
			if (bit > 0x80) {
				bit = 0x01;
				ptr_byte++;
			}
		}
	}
}

//---------------------------------------------------------------------------
// Function   : UpdateBext
// Description: bext チャンク更新
// ---
//	bext	:bext構造体へのアドレス
//	inFmt	:入力ファイル音声形式情報
//	inFmt	:入力ファイル音声形式情報
//	outFmt	:出力ファイル音声形式情報
//	param	:パラメーター構造体
//	bwf_size:bextのバイト数
//
int UpdateBext(BROADCAST_EXT *bext,SOUNDFMT *inFmt,SOUNDFMT *outFmt,PARAM_INFO *param,long bwf_size)
{
	__int64 timeRef;
	double d_timeRef;
	char strbuf[512];
	strbuf[0] = '\0';

	// version
	bext->version = 1;

	// TimeRef
	timeRef = bext->timeReferenceHigh;timeRef <<= 32;
	timeRef |= bext->timeReferenceLow;
	if (timeRef != 0) {
		d_timeRef = (double)timeRef;
		d_timeRef /= inFmt->sample;
		d_timeRef *= outFmt->sample;
		timeRef = (__int64)d_timeRef;
		bext->timeReferenceHigh = (DWORD)(timeRef >> 32);
		bext->timeReferenceLow	= (DWORD)timeRef;
	}
	// Coding History
	if (strlen(bext->codingHistory) == 0) {
		sprintf(strbuf,"A=PCM,F=%d,W=%d,M=%s,T=original,\r\n",inFmt->sample,inFmt->bitsPerSample,inFmt->channel == 1 ? "mono" : "stereo");
		strcat(bext->codingHistory,strbuf);
	}
	sprintf(strbuf,"A=PCM,F=%d,W=%d,M=%s,T=Upconv0.8.x;",outFmt->sample,outFmt->bitsPerSample,outFmt->channel == 1 ? "mono" : "stereo");
	strcat(bext->codingHistory,strbuf);
	if (param->abe) {
		strcpy(strbuf,"ABE;");
		strcat(bext->codingHistory,strbuf);
	}
	sprintf(strbuf,"HFA%d",param->hfa);
	strcat(bext->codingHistory,strbuf);
	strcat(bext->codingHistory,",\r\n");
	return 0;
}

//---------------------------------------------------------------------------
// Function   : dsf_main
// Description: 引数を処理し変換関数を呼び出す
//
//
int dsf_main(int argc, char *argv[])
{
	PARAM_INFO2 param;
	FILE *fp;
	char workpath[_MAX_PATH];
	char fracpath[_MAX_PATH];
	char infile[_MAX_PATH];
	char outfile[_MAX_PATH];
	char tmppath[_MAX_PATH];
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];
	char *p1,*p2;
	char param1[512];
	int retCode;
	long temp;

	do {
		infile[0] = '\0';
		outfile[0] = '\0';
		retCode = STATUS_SUCCESS;
		memset(&param,0,sizeof (PARAM_INFO2));

		if (argc < 3) {
			break;
		}

		strcpy(infile,argv[2]);
		strcpy(outfile,argv[3]);

		_splitpath(argv[3],drive,dir,fname,ext);
		_makepath(tmppath,drive,dir,fname,"param");
		// ファイルオープン
		fp = fopen(tmppath,"r");
		if (fp == NULL) {
			retCode = STATUS_FILE_READ_ERR;
			break;
		}

		if (fgets(workpath,_MAX_PATH,fp) == NULL) {
			retCode = STATUS_FILE_READ_ERR;
			break;
		}

		if (fgets(fracpath,_MAX_PATH,fp) == NULL) {
			retCode = STATUS_FILE_READ_ERR;
			break;
		}

		if (fgets(param1,512,fp) == NULL) {
			retCode = STATUS_FILE_READ_ERR;
			break;
		}
		p1 = strrchr(workpath,'\n');if (p1 != NULL) *p1 = '\0';
		if (strlen(workpath) >= 2 && workpath[strlen(workpath) - 1] != '\\') strcat(workpath,"\\");

		param.workpath = workpath;
		fclose(fp);
		p1 = param1;
		p2 = strchr(p1,(int)' ');

		for (;p1 != NULL;) {
			if (p2 != NULL) {
				*p2 = '\0';
			}

			if (p2 == NULL) {
				break;
			}
			p1 = p2 + 1;
			p2 = strchr(p1,(int)' ');
		}


	} while (0);

	if (param.err) {
		_splitpath(argv[3],drive,dir,fname,ext);
		_makepath(tmppath,drive,dir,fname,"err");
		fp = fopen(tmppath,"w");
		if (fp) {
			switch (param.err) {
				case STATUS_FILE_READ_ERR:
					fprintf(fp,"dsf2raw:[%04d] File read error.\n",param.errLine);
					break;
				case STATUS_FILE_WRITE_ERR:
					fprintf(fp,"dsf2raw:[%04d] File write error.\n",param.errLine);
					break;
				case STATUS_MEM_ALLOC_ERR:
					fprintf(fp,"dsf2raw:[%04d] Memory Allocation error.\n",param.errLine);
					break;
				default:
					fprintf(fp,"raw2wav:Other error.\n");
			}
			fclose(fp);
		}
	}
	return 0;
}

//---------------------------------------------------------------------------
// Function   : dsf_encode
// Description: DSF エンコード処理
// ---
// WAV ファイルを DSF ファイルへエンコードする
//
void dsf_encode(char *in_file,char *out_file,PARAM_INFO2 *param)
{
	// 未サポート
	param->err = 1;param->errLine = __LINE__;
}
//---------------------------------------------------------------------------
// Function   : dsf_decode
// Description: DSF デコード処理
// ---
// DSF ファイルを WAV ファイルへエンコードする
//
void dsf_decode(char *in_file,char *out_file,PARAM_INFO2 *param)
{
	char tmppath[_MAX_PATH];
	char workdrive[_MAX_DRIVE];
	char workdir[_MAX_DIR];
	char workfname[_MAX_FNAME];
	char workext[_MAX_EXT];
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];
	char tmpfile[_MAX_PATH];
	char work[128];
	FIO fp_r,fp_w1,fp_w2;
	FILE *fp;
	FILE *fp_files;
	DSF dsf;
	DSF_FMT dsf_fmt;
	DSF_DATA dsf_data;
	fio_size rs;
	fio_size ws;
	fio_size file_size;
	UI64 remain;
	char dbs[100];
	DWORD inSample,outSample;
	double dd;
	int lr_flag;
	int i,j,n;
	SSIZE avg;

#ifdef _OPENMP
	int nCpu;
	nCpu = param->thread;
	omp_set_num_threads(nCpu);
#endif

	do {
		memset(&fp_r,0,sizeof (FIO));
		memset(&fp_w1,0,sizeof (FIO));
		memset(&fp_w2,0,sizeof (FIO));
		fp_files = NULL;
		
		fio_open(&fp_r,in_file,FIO_MODE_R);
		if (fp_r.error) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		rs = fio_read(&dsf,sizeof (DSF),1,&fp_r);
		if (fp_r.error || rs != 1) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		
		if (memcmp(dsf.id,"DSD ",4)) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (dsf.chunk_size < 28) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (dsf.file_size < (28 + 52 + 12)) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		fio_seek(&fp_r,dsf.chunk_size,SEEK_SET);
		rs = fio_read(&dsf_fmt,sizeof (DSF_FMT),1,&fp_r);
		if (fp_r.error || rs != 1) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		
		if (memcmp(dsf_fmt.id,"fmt ",4)) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (dsf_fmt.chunk_size < 52) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (dsf_fmt.fmt_version != 1) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (dsf_fmt.fmt_id != 0) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (!(dsf_fmt.channel_type == 1 || dsf_fmt.channel_type == 2)) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (!(dsf_fmt.channel_count == 1 || dsf_fmt.channel_count == 2)) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (!(dsf_fmt.sampling == 2822400 || dsf_fmt.sampling == 2822400*2 || dsf_fmt.sampling == 2822400*4)) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (dsf_fmt.sample_bit_count != 1) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (dsf_fmt.sample_count < 1) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		if (dsf_fmt.block_size != 4096) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		fio_seek(&fp_r,dsf.chunk_size + dsf_fmt.chunk_size,SEEK_SET);
		rs = fio_read(&dsf_data,sizeof (DSF_DATA),1,&fp_r);
		if (fp_r.error || rs != 1) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		if (dsf_data.chunk_size <= 12) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		param->channel = dsf_fmt.channel_count;
		param->data_offset = fio_tell(&fp_r);

		sprintf(dbs,"DATA SIZE:%lld",dsf_data.chunk_size - 12);
		PRINT_LOG(dbs);
		sprintf(dbs,"sample count:%lld",dsf_fmt.sample_count);
		PRINT_LOG(dbs);

		_splitpath(out_file,drive,dir,fname,ext);
		_makepath(tmpfile,drive,dir,fname,"files");
		// 出力ファイルオープン
		fp_files = fopen(tmpfile,"a");
		if (fp_files == NULL) {
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}

		// 出力ファイル(Left)
		PRINT_LOG(out_file);
		_splitpath(out_file,drive,dir,fname,ext);
		_splitpath(param->workpath,workdrive,workdir,workfname,workext);
		PRINT_LOG(workdir);
		if (strlen(param->workpath) >= 3) {
			strcpy(workfname,fname);
		} else {
			strcpy(workdrive,drive);
			strcpy(workdir,dir);
			strcpy(workfname,fname);
		}

		_makepath(tmpfile,workdrive,workdir,workfname,"r1.tmp");
		PRINT_LOG(tmpfile);
		fprintf(fp_files,"%s\n",tmpfile);
		fio_open(&fp_w1,tmpfile,FIO_MODE_W);
		if (fp_w1.error) {
			PRINT_LOG("");
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}
		if (param->channel == 1) {
			fio_set_maxsize(&fp_w1,dsf_data.chunk_size - 12);
			fio_set_memory_limit(&fp_w1,20,param->fio);
		}
		if (param->channel == 2) {
			// 出力ファイル(Right)
			_makepath(tmpfile,workdrive,workdir,workfname,"r2.tmp");
			PRINT_LOG(tmpfile);
			fprintf(fp_files,"%s\n",tmpfile);
			fio_open(&fp_w2,tmpfile,FIO_MODE_W);
			if (fp_w2.error) {
				PRINT_LOG("");
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			fio_set_maxsize(&fp_w1,(dsf_data.chunk_size - 12) / 2);
			fio_set_maxsize(&fp_w2,(dsf_data.chunk_size - 12) / 2);
			fio_set_memory_limit(&fp_w1,20,param->fio);
			fio_set_memory_limit(&fp_w2,20,param->fio);
		}

		// 4096 ごとにインターリーブされているブロックをtmpファイルへ出力する(ステレオの場合)
		deinterleave(dsf_data.chunk_size - 12,&fp_r,&fp_w1,&fp_w2,param);
		if (param->err) {
			PRINT_LOG("");
			break;
		}
		fio_close(&fp_r);
		// 読み込みで開く
		fio_setmode_r(&fp_w1,&fp_r,NULL);
		if (fp_w1.error) {
			PRINT_LOG("");
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}
		if (fp_r.error) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		// 出力ファイル(Left)
		_makepath(tmpfile,workdrive,workdir,workfname,"r1");
		PRINT_LOG(tmpfile);
		fprintf(fp_files,"%s\n",tmpfile);
		fio_open(&fp_w1,tmpfile,FIO_MODE_W);
		if (fp_w1.error) {
			PRINT_LOG("");
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}

		param->inSampleR  = dsf_fmt.sampling;
//		param->outSampleR  = dsf_fmt.sampling;
		param->outSampleR = 192000;
		if (dsf_fmt.sampling == 2822400) {
			param->outSampleR = 192000;
			if (param->hfc == -1) param->hfc = 16000;
		} else if (dsf_fmt.sampling == 2822400*2) {
			param->outSampleR = 192000 * 2;
			if (param->hfc == -1) param->hfc = 22000;
		} else if (dsf_fmt.sampling == 2822400*4) {
			param->outSampleR = 192000 * 4;
			if (param->hfc == -1) param->hfc = 28000;
		}
		remain = dsf_fmt.sample_count;
		remain *= param->outSampleR;
		remain /= dsf_fmt.sampling;
		remain *= sizeof (SSIZE);

		// 出力ファイルサイズの最大値を指定
		fio_set_maxsize(&fp_w1,remain);
		fio_set_memory_limit(&fp_w1,20,param->fio);

		remain = dsf_fmt.sample_count;
		remain *= param->outSampleR;
		remain /= dsf_fmt.sampling;
		sprintf(dbs,"remain:%lld",remain);
		PRINT_LOG(dbs);
		param->n_sample = dsf_fmt.sample_count;

		// Left
		fftFilter(0,dsf_fmt.sample_count,remain,&fp_r,&fp_w1,param);
		if (param->err) {
			PRINT_LOG("");
			break;
		}
		fio_close(&fp_r);
		fio_close(&fp_w1);
		if (param->channel == 2) {
			// 読み込みで開く
			fio_setmode_r(&fp_w2,&fp_r,NULL);
			if (fp_w2.error) {
				PRINT_LOG("");
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			if (fp_r.error) {
				PRINT_LOG("");
				param->err = STATUS_FILE_READ_ERR;
				break;
			}

			// 出力ファイル(Right)
			_makepath(tmpfile,workdrive,workdir,workfname,"r2");
			PRINT_LOG(tmpfile);
			fprintf(fp_files,"%s\n",tmpfile);
			fio_open(&fp_w2,tmpfile,FIO_MODE_W);
			if (fp_w2.error) {
				PRINT_LOG("");
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			remain = dsf_fmt.sample_count;
			remain *= param->outSampleR;
			remain /= dsf_fmt.sampling;
			remain *= sizeof (SSIZE);

			// 出力ファイルサイズの最大値を指定
			fio_set_maxsize(&fp_w2,remain);
			fio_set_memory_limit(&fp_w2,20,param->fio);

			remain = dsf_fmt.sample_count;
			remain *= param->outSampleR;
			remain /= dsf_fmt.sampling;

			// Right
			fftFilter(1,dsf_fmt.sample_count,remain,&fp_r,&fp_w2,param);
			if (param->err) {
				PRINT_LOG("");
				break;
			}
			fio_close(&fp_r);
			fio_close(&fp_w2);
		}

		// 出力ファイル名の作成
		_splitpath(out_file,drive,dir,fname,ext);
		_makepath(tmpfile,drive,dir,fname,"param");
		// 出力ファイルオープン
		fp = fopen(tmpfile,"a");
		if (fp == NULL) {
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}
		sprintf(work," -is:%d -iw:%d -ch:%d",param->outSampleR,24,param->channel);
		strcat(param->argv4,work);
		dd = param->l_min;
		if (dd < 0) {
			dd *= -1;
		}
		if (dd < param->l_max) {
			dd = param->l_max;
		}
		
		dd /= (double)0x7FFFFFFFFFFFFF;

		if (param->l_cntLv > 0) {
			param->l_tmpLv /= param->l_cntLv;
		}

		if (param->l_maxLv > 0) {
			avg = (param->l_maxLv + param->l_tmpLv) / 2;
		} else {
			avg = param->l_tmpLv;
		}
		avg <<= 40;

		fprintf(fp,"r1=%.10lf,%llx\n",dd,avg);
//		fprintf(fp,"r1=%lf\n",1.00);
		sprintf(dbs,"r1=%lf",dd);
		PRINT_LOG(dbs);

		if (param->channel == 1) {
			fprintf(fp,"r2=%.10lf,%llx\n",1.0,0);
		} else {
			dd = param->r_min;
			if (dd < 0) {
				dd *= -1;
			}
			if (dd < param->r_max) {
				dd = param->r_max;
			}
			dd /= (double)0x7FFFFFFFFFFFFF;

			if (param->r_cntLv > 0) {
				param->r_tmpLv /= param->r_cntLv;
			}

			if (param->r_maxLv > 0) {
				avg = (param->r_maxLv + param->r_tmpLv) / 2;
			} else {
				avg = param->r_tmpLv;
			}
			avg <<= 40;

			fprintf(fp,"r2=%.10lf,%llx\n",dd,avg);
//			fprintf(fp,"r2=%lf\n",1.00);
			sprintf(dbs,"r2=%lf",dd);
			PRINT_LOG(dbs);
		}
		fprintf(fp,"r3=%lf,%llx\n",1.0,0);
		fprintf(fp,"r4=%lf,%llx\n",1.0,0);
		fprintf(fp,"r5=%lf,%llx\n",1.0,0);
		fprintf(fp,"r6=%lf,%llx\n",1.0,0);
		fclose(fp);
		fclose(fp_files);
	} while (0);
}
//---------------------------------------------------------------------------
// Function   : deinterleave
// Description: 4096ごとに記録されているチャンネルごとの1ブロックを結合する
// ---
//	fp_r		:入力ファイル用構造体
//	fp_w1		:出力ファイル用構造体(Left)
//	fp_w2		:出力ファイル用構造体(Right)
//	param		:変換パラメータ
//
void deinterleave(UI64 inByte,FIO *fp_r,FIO *fp_w1,FIO *fp_w2,PARAM_INFO2 *param)
{
	unsigned char *bit_buffer;
	fio_size rs,ws;
	int lr;

	bit_buffer	 = (unsigned char *)malloc(4096);
	if (bit_buffer == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
		return;
	}

	fio_seek(fp_r,param->data_offset,SEEK_SET);
	if (fp_r->error) {
		param->err = STATUS_FILE_READ_ERR;param->errLine = __LINE__;
		return;
	}

	lr = 0;
	while (inByte > 0) {
		rs = fio_read(bit_buffer,1,4096,fp_r);
		if (fp_r->error || rs != 4096) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;param->errLine = __LINE__;
			break;
		}
		if (lr == 0) {
			ws = fio_write(bit_buffer,1,4096,fp_w1);
			if (fp_w1->error || ws != rs) {
				PRINT_LOG("");
				param->err = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
				break;
			}
		} else {
			ws = fio_write(bit_buffer,1,4096,fp_w2);
			if (fp_w2->error || ws != rs) {
				PRINT_LOG("");
				param->err = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
				break;
			}
		}
		lr++;
		if (lr == 2) {
			lr = 0;
		}
		if (param->channel == 1) {
			lr = 0;
		}
		inByte -= 4096;
	}
	free(bit_buffer);
}
//---------------------------------------------------------------------------
// Function   : fftFilter
// Description: FFT によるフィルタ処理
// ---
//	lr			:Left/Right フラグ
//	inSample	:入力データのサンプル数(ch毎)
//	outSample	:出力データのサンプル数(ch毎)
//	fp_r		:入力ファイル用構造体
//	fp_w		:出力ファイル用構造体
//	param		:変換パラメータ
//
static void fftFilter(int lr,SSIZE inSample,SSIZE outSample,FIO *fp_r,FIO *fp_w,PARAM_INFO2 *param)
{
	SSIZE *mem0,*mem1;
	SSIZE *mem2,*mem3,*mem4;
	SSIZE *mem5,*mem6,*mem7;
	SSIZE *mem8;
	char *inPtr;
	long wkMemSize;
	long inSampleR,outSampleR;
	long fftSizeIn,fftSizeOut,i,j,n;
	long cutOff;
	long hfc;
	long wr;
	long wkSampleR;
	long pwCnt;
	double persent,per;
	double nx;
	double *pwBase,basePw;
	SSIZE *pIn[6],*pOut[6];
	SSIZE startInSample,inSampLen,outSampLen,nSample;
	UI64 outRemain;
	fftw_complex *fftw_io[7];
	fftw_plan fftw_p[7],fftw_ip[7];
	SSIZE membyte;
	SSIZE delay;
	SSIZE sigma;
	SSIZE a,b;
	SSIZE min,max;
	char s[50];
	SSIZE maxLv,maxLv2;
	SSIZE ns;

	fftw_io[0] = NULL;
	fftw_io[1] = NULL;
	fftw_io[2] = NULL;
	fftw_io[3] = NULL;
	fftw_io[4] = NULL;
	fftw_io[5] = NULL;
	fftw_io[6] = NULL;

	fio_rewind(fp_r);

	inSampleR = param->inSampleR;
	outSampleR = param->outSampleR;

	fftSizeIn = inSampleR / 5;
	fftSizeOut = outSampleR / 5;

	wkMemSize = fftSizeIn;
	wkMemSize *= 2;

	if (lr == 0) {
		min = param->l_min;
		max = param->l_max;
		maxLv  = param->l_maxLv;
		maxLv2 = param->l_tmpLv;
		ns	   = param->l_cntLv;
	} else {
		min = param->r_min;
		max = param->r_max;
		maxLv  = param->r_maxLv;
		maxLv2 = param->r_tmpLv;
		ns	   = param->r_cntLv;
	}
	
	// 入力用(1)
	mem1 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem1 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 出力用(1)
	mem2 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem2 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 出力用(1)
	mem3 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem3 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 出力用(1)
	mem4 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem4 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
#if 0
	// 出力用(2)
	mem5 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem5 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 出力用(2)
	mem6 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem6 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 出力用(2)
	mem7 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem7 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 出力用(3)
	mem8 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem8 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
#endif
	// 1
	fftw_io[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_io[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 2
	fftw_io[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_io[1] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 3
	fftw_io[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_io[2] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 4
	fftw_io[3] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_io[3] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
#if 0
	// 5
	fftw_io[4] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_io[4] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 6
	fftw_io[5] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_io[5] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 7
	fftw_io[6] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 192000 * 2);
	if (fftw_io[6] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
#endif

	// 2822400 → 192000 用のプラン(1)
	fftw_p[0] = fftw_plan_dft_1d(fftSizeIn,fftw_io[0],fftw_io[0],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[0] = fftw_plan_dft_1d(fftSizeOut,fftw_io[0],fftw_io[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 2822400 → 192000 用のプラン(2)
	fftw_p[1] = fftw_plan_dft_1d(fftSizeIn,fftw_io[1],fftw_io[1],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[1] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[1] = fftw_plan_dft_1d(fftSizeOut,fftw_io[1],fftw_io[1],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[1] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 2822400 → 192000 用のプラン(3)
	fftw_p[2] = fftw_plan_dft_1d(fftSizeIn,fftw_io[2],fftw_io[2],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[2] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[2] = fftw_plan_dft_1d(fftSizeOut,fftw_io[2],fftw_io[2],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[2] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
#if 0
	// 2822400 → 192000 用のプラン(4)
	fftw_p[3] = fftw_plan_dft_1d(fftSizeIn,fftw_io[3],fftw_io[3],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[3] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[3] = fftw_plan_dft_1d(fftSizeOut,fftw_io[3],fftw_io[3],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[3] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 2822400 → 192000 用のプラン(5)
	fftw_p[4] = fftw_plan_dft_1d(fftSizeIn,fftw_io[4],fftw_io[4],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[4] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[4] = fftw_plan_dft_1d(fftSizeOut,fftw_io[4],fftw_io[4],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[4] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 2822400 → 192000 用のプラン(6)
	fftw_p[5] = fftw_plan_dft_1d(fftSizeIn,fftw_io[5],fftw_io[5],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[5] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[5] = fftw_plan_dft_1d(fftSizeOut,fftw_io[5],fftw_io[5],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[5] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 192000 → 192000 用のプラン
	fftw_p[6] = fftw_plan_dft_1d(fftSizeOut,fftw_io[6],fftw_io[6],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[6] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[6] = fftw_plan_dft_1d(fftSizeOut,fftw_io[6],fftw_io[6],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[6] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
#endif

	outRemain = outSample;

	delay = 0;
	per = -1;
	sprintf(s,"%lld,%lld",inSample,outSample);
	PRINT_LOG(s);
//	for (startInSample = ((fftSizeIn + (fftSizeIn / 2)) * -1);startInSample < inSample + ((fftSizeIn * 2) + fftSizeIn / 2);startInSample += fftSizeIn * 2) {
	for (startInSample = ((fftSizeIn + (fftSizeIn / 2)) * -1);startInSample < inSample + ((fftSizeIn) + fftSizeIn / 2);startInSample += fftSizeIn) {
		if (startInSample >= 0 && startInSample < inSample) {
			persent = ((double)startInSample / inSample);
			persent *= 100;
			if (persent != per) {
				Sleep(1);
				if (param->channel == 2) {
					persent /= 2;
					if (lr == 1) {
						persent += 50;
					}
				}
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;
		}

		// 読み込んだデータをmem1に入れる
		memset(mem1,0,wkMemSize * sizeof (SSIZE));
		onebit2nbit(startInSample,wkMemSize,mem1,fp_r,param);
		if (param->err) {
			break;
		}

		// 1
		memset(mem2,0,wkMemSize * sizeof (SSIZE));
		memset(mem3,0,wkMemSize * sizeof (SSIZE));
		memset(mem4,0,wkMemSize * sizeof (SSIZE));

#if 0
		// 2
		memset(mem5,0,wkMemSize * sizeof (SSIZE));
		memset(mem6,0,wkMemSize * sizeof (SSIZE));
		memset(mem7,0,wkMemSize * sizeof (SSIZE));
		memset(mem8,0,wkMemSize * sizeof (SSIZE));
#endif
		// DSD -> 192000 へデシメーション(1)
		param->lfc = 2;
		if (param->dsf_mode != 1) {
			param->hfc = 96000;
		} else {
			param->hfc = 24000;
		}
		param->src_flag = 1;

PRINT_LOG("1");
		// スレッド1,2,3の組とスレッド4,5,6の組でfftする
		pIn[0]	= &mem1[((fftSizeIn / 2) * 0)];
		pOut[0] = &mem2[((fftSizeOut / 2) * 0)];
		pIn[1]	= &mem1[((fftSizeIn / 2) * 1)];
		pOut[1] = &mem3[((fftSizeOut / 2) * 1)];
		pIn[2]	= &mem1[((fftSizeIn / 2) * 2)];
		pOut[2] = &mem4[((fftSizeOut / 2) * 2)];

#if 0
		pIn[3]	= &mem1[((fftSizeIn / 2) * 2)];
		pOut[3] = &mem5[((fftSizeOut / 2) * 2)];
		pIn[4]	= &mem1[((fftSizeIn / 2) * 3)];
		pOut[4] = &mem6[((fftSizeOut / 2) * 3)];
		pIn[5]	= &mem1[((fftSizeIn / 2) * 4)];
		pOut[5] = &mem7[((fftSizeOut / 2) * 4)];
#endif
		#pragma omp parallel
		{
			#pragma omp sections
			{
				#pragma omp section
				{
					// 1
					fftFilterSub(pIn[0],pOut[0],fftw_io[0],fftw_io[0],fftw_p[0],fftw_ip[0],param,0);
				}
				#pragma omp section
				{
					// 2
					fftFilterSub(pIn[1],pOut[1],fftw_io[1],fftw_io[1],fftw_p[1],fftw_ip[1],param,1);
				}
				#pragma omp section
				{
					// 3
					fftFilterSub(pIn[2],pOut[2],fftw_io[2],fftw_io[2],fftw_p[2],fftw_ip[2],param,2);
				}
#if 0
				#pragma omp section
				{
					// 4
					fftFilterSub(pIn[3],pOut[3],fftw_io[3],fftw_io[3],fftw_p[3],fftw_ip[3],param,0);
				}
				#pragma omp section
				{
					// 5
					fftFilterSub(pIn[4],pOut[4],fftw_io[4],fftw_io[4],fftw_p[4],fftw_ip[4],param,1);
				}
				#pragma omp section
				{
					// 6
					fftFilterSub(pIn[5],pOut[5],fftw_io[5],fftw_io[5],fftw_p[5],fftw_ip[5],param,2);
				}
#endif
			}
		}
		if (param->dsf_mode != 1) {
			memset(mem1,0,fftSizeOut * 2 * sizeof (SSIZE));
			#pragma omp parallel for
			for (i = fftSizeOut / 2;i < fftSizeOut + fftSizeOut / 2;i++) {
				mem1[i] = mem2[i] + mem3[i] + mem4[i];
			}
#if 0
			#pragma omp parallel for
			for (i = fftSizeOut + fftSizeOut / 2;i < fftSizeOut * 2 + fftSizeOut / 2;i++) {
				mem1[i] = mem5[i] + mem6[i] + mem7[i];
			}
#endif

			for (i = fftSizeOut / 2;i < fftSizeOut * 2 + fftSizeOut / 2;i++) {
#if 0
				if (min > mem8[i]) {
					min = mem8[i];
				}
				if (max < mem8[i]) {
					max = mem8[i];
				}
#else
				if (min > mem1[i]) {
					min = mem1[i];
				}
				if (max < mem1[i]) {
					max = mem1[i];
				}
#endif
//				if ((mem8[i] >> 40) > 0) {
				if ((mem1[i] >> 40) > 0) {
//					maxLv2 += (mem8[i] >> 40);
					maxLv2 += (mem1[i] >> 40);
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


		} else {
//			memset(mem8,0,fftSizeOut * 3 * sizeof (SSIZE));
			memset(mem1,0,fftSizeOut * 2 * sizeof (SSIZE));
			#pragma omp parallel for
			for (i = fftSizeOut / 2;i < fftSizeOut + fftSizeOut / 2;i++) {
//				mem8[i] = mem2[i] + mem3[i] + mem4[i];
				mem1[i] = mem2[i] + mem3[i] + mem4[i];
			}
#if 0
			#pragma omp parallel for
			for (i = fftSizeOut + fftSizeOut / 2;i < fftSizeOut * 2 + fftSizeOut / 2;i++) {
//				mem8[i] = mem5[i] + mem6[i] + mem7[i];
				mem1[i] = mem5[i] + mem6[i] + mem7[i];
			}
#endif
			#pragma omp parallel for
			for (i = fftSizeOut / 2;i < fftSizeOut * 2 + fftSizeOut / 2;i++) {
#if 0
				if (min > mem8[i]) {
					min = mem8[i];
				}
				if (max < mem8[i]) {
					max = mem8[i];
				}
#else
				if (min > mem1[i]) {
					min = mem1[i];
				}
				if (max < mem1[i]) {
					max = mem1[i];
				}
#endif
//				if ((mem8[i] >> 40) > 0) {
				if ((mem1[i] >> 40) > 0) {
//					maxLv2 += (mem8[i] >> 40);
					maxLv2 += (mem1[i] >> 40);
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
#if 0
		if (param->dsf_mode != 1) {
			// DSD 24kz 以上のデータカット
	PRINT_LOG("3");
			param->lfc = -1;
			param->hfc = 24000;
			param->src_flag = 2;

			memset(mem2,0,wkMemSize * sizeof (SSIZE));
			memset(mem3,0,wkMemSize * sizeof (SSIZE));
			memset(mem4,0,wkMemSize * sizeof (SSIZE));
			memset(mem5,0,wkMemSize * sizeof (SSIZE));
			memset(mem6,0,wkMemSize * sizeof (SSIZE));
			memset(mem7,0,wkMemSize * sizeof (SSIZE));

			pIn[0]	= &mem1[((fftSizeOut / 2) * 0)];
			pOut[0] = &mem2[((fftSizeOut / 2) * 0)];
			pIn[1]	= &mem1[((fftSizeOut / 2) * 1)];
			pOut[1] = &mem3[((fftSizeOut / 2) * 1)];
			pIn[2]	= &mem1[((fftSizeOut / 2) * 2)];
			pOut[2] = &mem4[((fftSizeOut / 2) * 2)];
			pIn[3]	= &mem1[((fftSizeOut / 2) * 2)];
			pOut[3] = &mem5[((fftSizeOut / 2) * 2)];
			pIn[4]	= &mem1[((fftSizeOut / 2) * 3)];
			pOut[4] = &mem6[((fftSizeOut / 2) * 3)];
			pIn[5]	= &mem1[((fftSizeOut / 2) * 4)];
			pOut[5] = &mem7[((fftSizeOut / 2) * 4)];

			// 1
			fftFilterSub(pIn[0],pOut[0],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,0);
			// 2
			fftFilterSub(pIn[1],pOut[1],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,1);
			// 3
			fftFilterSub(pIn[2],pOut[2],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,2);
			// 4
			fftFilterSub(pIn[3],pOut[3],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,0);
			// 5
			fftFilterSub(pIn[4],pOut[4],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,1);
			// 6
			fftFilterSub(pIn[5],pOut[5],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,2);

			memset(mem8,0,fftSizeOut * 3 * sizeof (SSIZE));
			#pragma omp parallel for
			for (i = fftSizeOut / 2;i < fftSizeOut + fftSizeOut / 2;i++) {
				mem8[i] = mem2[i] + mem3[i] + mem4[i];
			}
			#pragma omp parallel for
			for (i = fftSizeOut + fftSizeOut / 2;i < fftSizeOut * 2 + fftSizeOut / 2;i++) {
				mem8[i] = mem5[i] + mem6[i] + mem7[i];
			}
			#pragma omp parallel for
			for (i = fftSizeOut / 2;i < fftSizeOut * 2 + fftSizeOut / 2;i++) {
				if (min > mem8[i]) {
					min = mem8[i];
				}
				if (max < mem8[i]) {
					max = mem8[i];
				}
				if ((mem8[i] >> 40) > 0) {
					maxLv2 += (mem8[i] >> 40);
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

			// 24kHz 以下のデータカット
			param->lfc = 24000;
			param->hfc = -1;
			param->src_flag = 2;

			memset(mem2,0,wkMemSize * sizeof (SSIZE));
			memset(mem3,0,wkMemSize * sizeof (SSIZE));
			memset(mem4,0,wkMemSize * sizeof (SSIZE));
			memset(mem5,0,wkMemSize * sizeof (SSIZE));
			memset(mem6,0,wkMemSize * sizeof (SSIZE));
			memset(mem7,0,wkMemSize * sizeof (SSIZE));

			pIn[0]	= &mem1[((fftSizeOut / 2) * 0)];
			pOut[0] = &mem2[((fftSizeOut / 2) * 0)];
			pIn[1]	= &mem1[((fftSizeOut / 2) * 1)];
			pOut[1] = &mem3[((fftSizeOut / 2) * 1)];
			pIn[2]	= &mem1[((fftSizeOut / 2) * 2)];
			pOut[2] = &mem4[((fftSizeOut / 2) * 2)];
			pIn[3]	= &mem1[((fftSizeOut / 2) * 2)];
			pOut[3] = &mem5[((fftSizeOut / 2) * 2)];
			pIn[4]	= &mem1[((fftSizeOut / 2) * 3)];
			pOut[4] = &mem6[((fftSizeOut / 2) * 3)];
			pIn[5]	= &mem1[((fftSizeOut / 2) * 4)];
			pOut[5] = &mem7[((fftSizeOut / 2) * 4)];

			// 1
			fftFilterSub(pIn[0],pOut[0],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,0);
			// 2
			fftFilterSub(pIn[1],pOut[1],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,1);
			// 3
			fftFilterSub(pIn[2],pOut[2],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,2);
			// 4
			fftFilterSub(pIn[3],pOut[3],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,0);
			// 5
			fftFilterSub(pIn[4],pOut[4],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,1);
			// 6
			fftFilterSub(pIn[5],pOut[5],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,2);


			memset(mem1,0,fftSizeOut * 2 * sizeof (SSIZE));
			#pragma omp parallel for
			for (i = fftSizeOut / 2;i < fftSizeOut + fftSizeOut / 2;i++) {
				mem1[i] = mem2[i] + mem3[i] + mem4[i];
			}
			#pragma omp parallel for
			for (i = fftSizeOut + fftSizeOut / 2;i < fftSizeOut * 2 + fftSizeOut / 2;i++) {
				mem1[i] = mem5[i] + mem6[i] + mem7[i];
			}

			if (param->dsf_mode == 2) {
				ana_abe(0,fftSizeOut * 3,mem1,mem5,param);
				#pragma omp parallel for
				for (i = 0;i < fftSizeOut * 3;i++) {
					mem1[i] = mem5[i];
				}
			}
			// 24kHz 以下のデータカット
			param->lfc = 24000;
			param->hfc = -1;
			param->src_flag = 2;

			memset(mem2,0,wkMemSize * sizeof (SSIZE));
			memset(mem3,0,wkMemSize * sizeof (SSIZE));
			memset(mem4,0,wkMemSize * sizeof (SSIZE));
			memset(mem5,0,wkMemSize * sizeof (SSIZE));
			memset(mem6,0,wkMemSize * sizeof (SSIZE));
			memset(mem7,0,wkMemSize * sizeof (SSIZE));

			pIn[0]	= &mem1[((fftSizeOut / 2) * 0)];
			pOut[0] = &mem2[((fftSizeOut / 2) * 0)];
			pIn[1]	= &mem1[((fftSizeOut / 2) * 1)];
			pOut[1] = &mem3[((fftSizeOut / 2) * 1)];
			pIn[2]	= &mem1[((fftSizeOut / 2) * 2)];
			pOut[2] = &mem4[((fftSizeOut / 2) * 2)];
			pIn[3]	= &mem1[((fftSizeOut / 2) * 2)];
			pOut[3] = &mem5[((fftSizeOut / 2) * 2)];
			pIn[4]	= &mem1[((fftSizeOut / 2) * 3)];
			pOut[4] = &mem6[((fftSizeOut / 2) * 3)];
			pIn[5]	= &mem1[((fftSizeOut / 2) * 4)];
			pOut[5] = &mem7[((fftSizeOut / 2) * 4)];

			// 1
			fftFilterSub(pIn[0],pOut[0],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,0);
			// 2
			fftFilterSub(pIn[1],pOut[1],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,1);
			// 3
			fftFilterSub(pIn[2],pOut[2],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,2);
			// 4
			fftFilterSub(pIn[3],pOut[3],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,0);
			// 5
			fftFilterSub(pIn[4],pOut[4],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,1);
			// 6
			fftFilterSub(pIn[5],pOut[5],fftw_io[6],fftw_io[6],fftw_p[6],fftw_ip[6],param,2);


			#pragma omp parallel for
			for (i = fftSizeOut / 2;i < fftSizeOut + fftSizeOut / 2;i++) {
				mem8[i] += mem2[i] + mem3[i] + mem4[i];
			}
			#pragma omp parallel for
			for (i = fftSizeOut + fftSizeOut / 2;i < fftSizeOut * 2 + fftSizeOut / 2;i++) {
				mem8[i] += mem5[i] + mem6[i] + mem7[i];
			}
		}
#endif
		if (startInSample + fftSizeIn / 2 >= 0) {
			if (outRemain >= fftSizeOut) {
PRINT_LOG("6");
//				wr = fio_write(mem8 + (fftSizeOut / 2),sizeof (SSIZE),fftSizeOut * 2,fp_w);
				wr = fio_write(mem1 + (fftSizeOut / 2),sizeof (SSIZE),fftSizeOut,fp_w);
				if (wr != fftSizeOut) {
					char s[100];
					sprintf(s,"%ld:fio_write(%ld,%ld)",wr,sizeof (SSIZE),fftSizeOut);
					param->err = STATUS_FILE_WRITE_ERR;
					PRINT_LOG(s);
					return;
				}
			} else {
//				wr = fio_write(mem8 + (fftSizeOut / 2),sizeof (SSIZE),outRemain,fp_w);
				wr = fio_write(mem1 + (fftSizeOut / 2),sizeof (SSIZE),outRemain,fp_w);
				if (wr != outRemain) {
					param->err = STATUS_FILE_WRITE_ERR;
					PRINT_LOG("ERROR");
					return;
				}
			}
			if (outRemain >= fftSizeOut) {
				outRemain -= fftSizeOut;
			} else {
				break;
			}
		}
PRINT_LOG("7");
	}
	PRINT_LOG("end");

	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	al_free(mem4);
#if 0
	al_free(mem5);
	al_free(mem6);
	al_free(mem7);
	al_free(mem8);
#endif
	// 1
	fftw_destroy_plan(fftw_p[0]);
	fftw_destroy_plan(fftw_ip[0]);
	fftw_destroy_plan(fftw_p[1]);
	fftw_destroy_plan(fftw_ip[1]);
	fftw_destroy_plan(fftw_p[2]);
	fftw_destroy_plan(fftw_ip[2]);

#if 0
	fftw_destroy_plan(fftw_p[3]);
	fftw_destroy_plan(fftw_ip[3]);
	fftw_destroy_plan(fftw_p[4]);
	fftw_destroy_plan(fftw_ip[4]);
	fftw_destroy_plan(fftw_p[5]);
	fftw_destroy_plan(fftw_ip[5]);
	fftw_destroy_plan(fftw_p[6]);
	fftw_destroy_plan(fftw_ip[6]);
#endif
	fftw_free(fftw_io[0]);
	fftw_free(fftw_io[1]);
	fftw_free(fftw_io[2]);
#if 0
	fftw_free(fftw_io[3]);
	fftw_free(fftw_io[4]);
	fftw_free(fftw_io[5]);
	fftw_free(fftw_io[6]);
#endif
	if (lr == 0) {
		param->l_min = min;
		param->l_max = max;
		param->l_maxLv = maxLv;
		param->l_tmpLv = maxLv2;
		param->l_cntLv = ns;
	} else {
		param->r_min = min;
		param->r_max = max;
		param->r_maxLv = maxLv;
		param->r_tmpLv = maxLv2;
		param->r_cntLv = ns;
	}
}
//---------------------------------------------------------------------------
// Function   : fftFilterSub
// Description: FFT によるフィルタ処理(サブ関数)
// ---
//	param		:変換パラメータ
//
static void fftFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,PARAM_INFO2 *param,int id)
{
	long inSampleR,outSampleR;
	long wkSampleR;
	long fftSizeIn,fftSizeOut,i,j,n;
	long cutOff;
	long hfc,lfc;
	double nx;
	double p;
	long validIndex;
	long index1k,index15k,index19k,width1k,width_l,width_h,index_l,index_h,index_b;
	long h;
	double pw;
	double rate;
	double wid;
PRINT_LOG("Start");
	inSampleR = param->inSampleR;
	outSampleR = param->outSampleR;
	if (param->src_flag == 2) {
		inSampleR = param->outSampleR;
	}

	fftSizeIn = inSampleR / 5;
	fftSizeOut = outSampleR / 5;

	for (i = 0;i < fftSizeOut;i++) {
		fftw_out[i][0] = 0;
		fftw_out[i][1] = 0;
	}

	// FFT 初期設定
	copyToFFTW(fftw_in,pIn,fftSizeIn);

	windowFFTW(fftw_in,fftSizeIn);
	
	// FFT
	fftw_execute(fftw_p);

	// 高域削除
	if (inSampleR <= outSampleR) {
		wkSampleR = inSampleR;
	} else {
		wkSampleR = outSampleR;
	}
	hfc = wkSampleR / 2;
	
	if (param->hfc != -1) {
		if (hfc > param->hfc) {
			hfc = param->hfc;
		}
	}
	cutOff = ((double)fftSizeOut / outSampleR) * hfc;
	cutFFTW(fftw_out,cutOff,fftSizeOut);

	if (param->lfc != -1) {
		cutOff = ((double)fftSizeOut / outSampleR) * param->lfc;
		for (i = 1;i < cutOff;i++) {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
	}

	// 半分のデータを復元
	for (i = 1;i < fftSizeOut / 2;i++) {
		fftw_out[fftSizeOut - i][0] = fftw_out[i][0];
		fftw_out[fftSizeOut - i][1] = fftw_out[i][1] * -1;
	}

	fftw_out[0][0] = 0;
	fftw_out[0][1] = 0;
	fftw_out[1][0] = 0;
	fftw_out[1][1] = 0;

	// invert FFT
	fftw_execute(fftw_ip);

	// 出力
	for (i = 0;i < fftSizeOut;i++) {
		pOut[i] = (SSIZE)(fftw_in[i][0] / fftSizeOut);
	}
PRINT_LOG("End");
}
//---------------------------------------------------------------------------
// Function   : onebit2nbit
// Description: 1bit データを 64bit データへ変換して読み込む
// ---
//	offset		:入力オフセット(サンプル毎)
//	n			:サンプル数
//	buffer		:データバッファ
//	fp_r		:入力FIO
//	param		:変換パラメータ
//
void onebit2nbit(SSIZE offset,SSIZE sample,SSIZE *buffer,FIO *fp_r,PARAM_INFO2 *param)
{
	unsigned char *bit_buffer;
	int *delay_buffer;
	unsigned char mask;
	double nbit_r;
	SSIZE nbit_data;
	SSIZE seek_ptr;
	int i,j,n,bit_count_p,bit_count_m;
	int ptr_in_byte,ptr_out_byte;
	int ptr_delay;
	int delay;
	int n_delay;
	int zero_count;
	fio_size rs;
	char sss[50];
	
	PRINT_LOG("Start");

	memset(buffer,0,sample * sizeof (SSIZE));
	if (offset < 0) {
		if ((offset * -1) >= sample) {
			// ファイルから読み込むデータがないのでリターンする
			PRINT_LOG("No Data");
			return;
		}
		buffer += (offset * -1);
		sample -= (offset * -1);
		offset	= 0;
	}

	if (offset + sample > param->n_sample) {
		sample = param->n_sample - offset;
		sprintf(sss,"offset:%lld,n:%lld",offset,sample);
		PRINT_LOG(sss);
	}

	if (sample <= 0) {
		return;
	}

	delay = 24;

	bit_buffer	 = (unsigned char *)malloc(4096);
	delay_buffer = (int *)malloc(delay * 2 * sizeof (int));
	if (bit_buffer == NULL || delay_buffer == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	ptr_delay = delay;
	
	// delay_buffer には3値を使う(1,-1,0)
	for (i = 0;i < delay * 2;i++) {
		delay_buffer[i] = 0;	// 無音
	}

	seek_ptr = offset - delay;
	n_delay = delay;
	if ((seek_ptr + delay) >= 0) {
		// delay バッファに以前の値をセットする
		if (seek_ptr < 0) {
			n_delay += seek_ptr;
			seek_ptr = 0;
		}

		ptr_out_byte = 0;
		while (n_delay > 0) {
			fio_seek(fp_r,(seek_ptr / 8) + param->data_offset,SEEK_SET);
			if (fp_r->error) {
				param->err = STATUS_FILE_READ_ERR;
				return;
			}
			rs = fio_read(bit_buffer,1,(n_delay + 7) / 8,fp_r);
			if (fp_r->error) {
				PRINT_LOG("");
				param->err = STATUS_FILE_READ_ERR;
				return;
			}
			if (rs == 0) {
				break;
			}
			ptr_in_byte = 0;
			mask = 0x01 << (seek_ptr % 8);
//			if (bit_buffer[ptr_in_byte] == 0x69 && bit_buffer[ptr_in_byte + 1] == 0x69 && bit_buffer[ptr_in_byte + 2] == 0x69) {
				for (i = 0;i < 8 * 3;i++) {
					delay_buffer[i] = 0;
				}
				n_delay = 0;
//			} else {
				do {
					if ((bit_buffer[ptr_in_byte] & mask)) {
						delay_buffer[ptr_out_byte] = 1;
					} else {
						delay_buffer[ptr_out_byte] = -1;
					}
					mask <<= 1;
					if (mask == 0x00) {
						mask = 0x01;
						ptr_in_byte++;
						rs--;
					}
					ptr_out_byte++;
					seek_ptr++;
					n_delay--;
				} while (rs && n_delay > 0);
//			}
		}
	}

	bit_count_p = bit_count_m = -1;
	ptr_out_byte = 0;
	zero_count = 0;

	while (sample > 0) {
		seek_ptr = offset;
		seek_ptr /= 8;			// bit
		fio_seek(fp_r,seek_ptr,SEEK_SET);
		if (fp_r->error) {
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		rs = fio_read(bit_buffer,1,4096,fp_r);
		if (fp_r->error) {
			PRINT_LOG("");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		seek_ptr = offset;
		seek_ptr /= 8;			// bit
		ptr_in_byte = 0;
		mask = 0x01 << (offset % 8);

		do {
			if (bit_count_p != -1 && bit_count_m != -1) {
				if (delay_buffer[ptr_delay - delay] == 1) {
					bit_count_p--;
				} else if (delay_buffer[ptr_delay - delay] == -1) {
					bit_count_m--;
				}
			}
#if 0
			if (zero_count == 0 && ptr_in_byte + 2 < rs && bit_buffer[ptr_in_byte + 0] == 0x69 && bit_buffer[ptr_in_byte + 1] == 0x69 && bit_buffer[ptr_in_byte + 2] == 0x69) {
				zero_count = 24;
				zero_count -= (offset % 8);
			} else if (zero_count == 0 && ptr_in_byte + 2 == rs && bit_buffer[ptr_in_byte + 0] == 0x69 && bit_buffer[ptr_in_byte + 1] == 0x69) {
				zero_count = 16;
				zero_count -= (offset % 8);
			} else if (zero_count == 0 && sample > 0 && ptr_in_byte + 1 == rs && bit_buffer[ptr_in_byte + 0] == 0x69) {
				zero_count = 8;
				zero_count -= (offset % 8);
			}
#endif
			if (zero_count == 0) {
				if ((bit_buffer[ptr_in_byte] & mask)) {
					delay_buffer[ptr_delay] = 1;
					delay_buffer[ptr_delay - delay] = 1;
					if (bit_count_p != -1) {
						bit_count_p++;
					}
				} else {
					delay_buffer[ptr_delay] = -1;
					delay_buffer[ptr_delay - delay] = -1;
					if (bit_count_m != -1) {
						bit_count_m++;
					}
				}
			} else {
				delay_buffer[ptr_delay] = 0;
				delay_buffer[ptr_delay - delay] = 0;
			}
			if (bit_count_p == -1 && bit_count_m == -1) {
				for (i = 0,bit_count_p = 0,bit_count_m = 0;i < delay;i++) {
					if (delay_buffer[ptr_delay - i] == 1) {
						bit_count_p++;
					} else if (delay_buffer[ptr_delay - i] == -1) {
						bit_count_m++;
					}
				}
			}
			if ((bit_count_p + bit_count_m) > 0) {
				nbit_r = (double)1.0 * ((double)bit_count_p / (bit_count_p + bit_count_m));
				if ((bit_count_p + bit_count_m) < delay) {
					nbit_r = (nbit_r / delay) * (bit_count_p + bit_count_m) + ((double)0.5 / delay) * (delay - (bit_count_p + bit_count_m));
				}
				if (param->inSampleR == 2822400) {
					nbit_data = (SSIZE)1 << 53;
					nbit_data = (SSIZE)((double)nbit_data * nbit_r);
					nbit_data -= (SSIZE)1 << (53 - 1);
				} else if (param->inSampleR == 2822400 * 2) {
					nbit_data = (SSIZE)1 << 51;
					nbit_data = (SSIZE)((double)nbit_data * nbit_r);
					nbit_data -= (SSIZE)1 << (51 - 1);
				} else {
					nbit_data = (SSIZE)1 << 49;
					nbit_data = (SSIZE)((double)nbit_data * nbit_r);
					nbit_data -= (SSIZE)1 << (49 - 1);
				}
			} else {
				if (param->inSampleR == 2822400) {
					nbit_data = (SSIZE)1 << 53;
					nbit_data = (SSIZE)((double)nbit_data * 0.5);
					nbit_data -= (SSIZE)1 << (53 - 1);
				} else if (param->inSampleR == 2822400 * 2) {
					nbit_data = (SSIZE)1 << 51;
					nbit_data = (SSIZE)((double)nbit_data * 0.5);
					nbit_data -= (SSIZE)1 << (51 - 1);
				} else {
					nbit_data = (SSIZE)1 << 49;
					nbit_data = (SSIZE)((double)nbit_data * 0.5);
					nbit_data -= (SSIZE)1 << (49 - 1);
				}
			}
			if (zero_count > 0) {
				zero_count--;
			}
			if (sample > 0) {
				buffer[ptr_out_byte] = nbit_data;
				ptr_out_byte ++;
				offset++;
				sample--;
			}
			ptr_delay++;
			if (ptr_delay >= delay*2) {
				ptr_delay = delay;
			}
			mask <<= 1;
			if (mask == 0) {
				mask = 0x01;
				ptr_in_byte++;
			}
		} while (ptr_in_byte < rs);
	}
	free(bit_buffer);
	free(delay_buffer);
	PRINT_LOG("End");
}
//---------------------------------------------------------------------------
// Function   : ana_abe
// Description: 高域ノイズ軽減処理
// ---
//	start		:開始位置
//	nSample		:サンプル数
//	i_buffer	:入力バッファ
//	o_buffer	:出力バッファ
//	param		:パラメーター構造体
//	
//
void ana_abe(SSIZE start,UI64 nSample,SSIZE *i_buffer,SSIZE *o_buffer,PARAM_INFO2 *param)
{
	int i,j,n;
	
	SSIZE min,max,dd;
PRINT_LOG("Start");
	if (start < 3) {
		for (i = start + 3,j = 0;j < nSample;i++,j++) {
			for (n = 0,min = i_buffer[i],max = i_buffer[i];n < 3;n++) {
				if (min > i_buffer[i - n]) {
					min = i_buffer[i - n];
				} else if (max < i_buffer[i - n]) {
					max = i_buffer[i - n];
				}
			}
			o_buffer[i] = (min + max) / 2;
		}
	} else {
		for (i = start,j = 0;j < nSample;i++,j++) {
			for (n = 0,min = i_buffer[i],max = i_buffer[i];n < 3;n++) {
				if (min > i_buffer[i - n]) {
					min = i_buffer[i - n];
				} else if (max < i_buffer[i - n]) {
					max = i_buffer[i - n];
				}
			}
			o_buffer[i] = (min + max) / 2;
		}
	}
	if (param->inSampleR == 2822400) {
		if (start < 10) {
			for (i = start + 10,j = 0;j < nSample;i++,j++) {
				for (n = 0,min = o_buffer[i],max = o_buffer[i];n < 10;n++) {
					if (min > o_buffer[i - n]) {
						min = o_buffer[i - n];
					} else if (max < o_buffer[i - n]) {
						max = o_buffer[i - n];
					}
				}
				o_buffer[i] = (min + max) / 2;
			}
		} else {
			for (i = start,j = 0;j < nSample;i++,j++) {
				for (n = 0,min = o_buffer[i],max = o_buffer[i];n < 10;n++) {
					if (min > o_buffer[i - n]) {
						min = o_buffer[i - n];
					} else if (max < o_buffer[i - n]) {
						max = o_buffer[i - n];
					}
				}
				o_buffer[i] = (min + max) / 2;
			}
		}

		for (i = start + 1,j = 0;j + 1 < nSample;i++,j++) {
			o_buffer[i] = (o_buffer[i - 1] + o_buffer[i + 1]) / 2;
		}
		for (i = start + 1,j = 0;j + 1 < nSample;i++,j++) {
			o_buffer[i] = (o_buffer[i - 1] + o_buffer[i + 1]) / 2;
		}
	}
#if 0
	for (i = start + 1,j = 0;j + 1 < nSample;i++,j++) {
		o_buffer[i] = (o_buffer[i - 1] + o_buffer[i + 1]) / 2;
	}
	for (i = start + 1,j = 0;j + 1 < nSample;i++,j++) {
		o_buffer[i] = (o_buffer[i - 1] + o_buffer[i + 1]) / 2;
	}
	for (i = start + 1,j = 0;j + 1 < nSample;i++,j++) {
		o_buffer[i] = (o_buffer[i - 1] + o_buffer[i + 1]) / 2;
	}
	for (i = start + 1,j = 0;j + 1 < nSample;i++,j++) {
		o_buffer[i] = (o_buffer[i - 1] + o_buffer[i + 1]) / 2;
	}
#endif
	PRINT_LOG("End");
}
#if 0
//---------------------------------------------------------------------------
// Function   : copyToFFTW
// Description: fftw用配列に値をコピーする
// ---
//	
//
void copyToFFTW(fftw_complex *fftw,SSIZE *buf,long size)
{
	long i;
	
	for (i = 0;i + 64 < size;i+=64) {
		fftw[i + 0][0] = buf[0];
		fftw[i + 0][1] = 0;
		fftw[i + 1][0] = buf[1];
		fftw[i + 1][1] = 0;
		fftw[i + 2][0] = buf[2];
		fftw[i + 2][1] = 0;
		fftw[i + 3][0] = buf[3];
		fftw[i + 3][1] = 0;
		fftw[i + 4][0] = buf[4];
		fftw[i + 4][1] = 0;
		fftw[i + 5][0] = buf[5];
		fftw[i + 5][1] = 0;
		fftw[i + 6][0] = buf[6];
		fftw[i + 6][1] = 0;
		fftw[i + 7][0] = buf[7];
		fftw[i + 7][1] = 0;
		fftw[i + 8][0] = buf[8];
		fftw[i + 8][1] = 0;
		fftw[i + 9][0] = buf[9];
		fftw[i + 9][1] = 0;
		fftw[i + 10][0] = buf[10];
		fftw[i + 10][1] = 0;
		fftw[i + 11][0] = buf[11];
		fftw[i + 11][1] = 0;
		fftw[i + 12][0] = buf[12];
		fftw[i + 12][1] = 0;
		fftw[i + 13][0] = buf[13];
		fftw[i + 13][1] = 0;
		fftw[i + 14][0] = buf[14];
		fftw[i + 14][1] = 0;
		fftw[i + 15][0] = buf[15];
		fftw[i + 15][1] = 0;
		fftw[i + 16][0] = buf[16];
		fftw[i + 16][1] = 0;
		fftw[i + 17][0] = buf[17];
		fftw[i + 17][1] = 0;
		fftw[i + 18][0] = buf[18];
		fftw[i + 18][1] = 0;
		fftw[i + 19][0] = buf[19];
		fftw[i + 19][1] = 0;
		fftw[i + 20][0] = buf[20];
		fftw[i + 20][1] = 0;
		fftw[i + 21][0] = buf[21];
		fftw[i + 21][1] = 0;
		fftw[i + 22][0] = buf[22];
		fftw[i + 22][1] = 0;
		fftw[i + 23][0] = buf[23];
		fftw[i + 23][1] = 0;
		fftw[i + 24][0] = buf[24];
		fftw[i + 24][1] = 0;
		fftw[i + 25][0] = buf[25];
		fftw[i + 25][1] = 0;
		fftw[i + 26][0] = buf[26];
		fftw[i + 26][1] = 0;
		fftw[i + 27][0] = buf[27];
		fftw[i + 27][1] = 0;
		fftw[i + 28][0] = buf[28];
		fftw[i + 28][1] = 0;
		fftw[i + 29][0] = buf[29];
		fftw[i + 29][1] = 0;
		fftw[i + 30][0] = buf[30];
		fftw[i + 30][1] = 0;
		fftw[i + 31][0] = buf[31];
		fftw[i + 31][1] = 0;
		fftw[i + 32][0] = buf[32];
		fftw[i + 32][1] = 0;
		fftw[i + 33][0] = buf[33];
		fftw[i + 33][1] = 0;
		fftw[i + 34][0] = buf[34];
		fftw[i + 34][1] = 0;
		fftw[i + 35][0] = buf[35];
		fftw[i + 35][1] = 0;
		fftw[i + 36][0] = buf[36];
		fftw[i + 36][1] = 0;
		fftw[i + 37][0] = buf[37];
		fftw[i + 37][1] = 0;
		fftw[i + 38][0] = buf[38];
		fftw[i + 38][1] = 0;
		fftw[i + 39][0] = buf[39];
		fftw[i + 39][1] = 0;
		fftw[i + 40][0] = buf[40];
		fftw[i + 40][1] = 0;
		fftw[i + 41][0] = buf[41];
		fftw[i + 41][1] = 0;
		fftw[i + 42][0] = buf[42];
		fftw[i + 42][1] = 0;
		fftw[i + 43][0] = buf[43];
		fftw[i + 43][1] = 0;
		fftw[i + 44][0] = buf[44];
		fftw[i + 44][1] = 0;
		fftw[i + 45][0] = buf[45];
		fftw[i + 45][1] = 0;
		fftw[i + 46][0] = buf[46];
		fftw[i + 46][1] = 0;
		fftw[i + 47][0] = buf[47];
		fftw[i + 47][1] = 0;
		fftw[i + 48][0] = buf[48];
		fftw[i + 48][1] = 0;
		fftw[i + 49][0] = buf[49];
		fftw[i + 49][1] = 0;
		fftw[i + 50][0] = buf[50];
		fftw[i + 50][1] = 0;
		fftw[i + 51][0] = buf[51];
		fftw[i + 51][1] = 0;
		fftw[i + 52][0] = buf[52];
		fftw[i + 52][1] = 0;
		fftw[i + 53][0] = buf[53];
		fftw[i + 53][1] = 0;
		fftw[i + 54][0] = buf[54];
		fftw[i + 54][1] = 0;
		fftw[i + 55][0] = buf[55];
		fftw[i + 55][1] = 0;
		fftw[i + 56][0] = buf[56];
		fftw[i + 56][1] = 0;
		fftw[i + 57][0] = buf[57];
		fftw[i + 57][1] = 0;
		fftw[i + 58][0] = buf[58];
		fftw[i + 58][1] = 0;
		fftw[i + 59][0] = buf[59];
		fftw[i + 59][1] = 0;
		fftw[i + 60][0] = buf[60];
		fftw[i + 60][1] = 0;
		fftw[i + 61][0] = buf[61];
		fftw[i + 61][1] = 0;
		fftw[i + 62][0] = buf[62];
		fftw[i + 62][1] = 0;
		fftw[i + 63][0] = buf[63];
		fftw[i + 63][1] = 0;
		buf += 64;
	}
	for (;i < size;i++) {
		fftw[i][0] = *buf++;
		fftw[i][1] = 0;
	}
}
//---------------------------------------------------------------------------
// Function   : windowFFTW
// Description: FFTW用Window関数
// ---
//	
//
void windowFFTW(fftw_complex *fftw,long size)
{
	long i,j;

	// ウインドウサイズ毎に定数化する
	switch (size) {
		case (4096 * 1):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 1) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 1));
			}
			#pragma omp parallel for
			for (i = ((4096 * 1) - 1) / 2;i < (4096 * 1);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 1));
			}
			break;
		case (4096 * 2):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 2) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 2));
			}
			#pragma omp parallel for
			for (i = ((4096 * 2) - 1) / 2;i < (4096 * 2);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 2));
			}
			break;
		case (4096 * 4):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 4) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 4));
			}
			#pragma omp parallel for
			for (i = ((4096 * 4) - 1) / 2;i < (4096 * 4);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 4));
			}
			break;
		case (4096 * 8):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 8) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 8));
			}
			#pragma omp parallel for
			for (i = ((4096 * 8) - 1) / 2;i < (4096 * 8);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 8));
			}
			break;
		case (4096 * 16):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 16) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 16));
			}
			#pragma omp parallel for
			for (i = ((4096 * 16) - 1) / 2;i < (4096 * 16);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 16));
			}
			break;
		case (32000 * 1):
			#pragma omp parallel for
			for (i = 0;i < ((32000 * 1) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(32000 * 1));
			}
			#pragma omp parallel for
			for (i = ((32000 * 1) - 1) / 2;i < (32000 * 1);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(32000 * 1));
			}
			break;
		case (44100 * 1):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 1) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 1));
			}
			#pragma omp parallel for
			for (i = ((44100 * 1) - 1) / 2;i < (44100 * 1);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 1));
			}
			break;
		case (48000 * 1):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 1) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 1));
			}
			#pragma omp parallel for
			for (i = ((48000 * 1) - 1) / 2;i < (48000 * 1);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 1));
			}
			break;
		case (44100 * 2):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 2) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 2));
			}
			#pragma omp parallel for
			for (i = ((44100 * 2) - 1) / 2;i < (44100 * 2);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 2));
			}
			break;
		case (48000 * 2):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 2) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 2));
			}
			#pragma omp parallel for
			for (i = ((48000 * 2) - 1) / 2;i < (48000 * 2);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 2));
			}
			break;
		case (44100 * 4):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 4) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 4));
			}
			#pragma omp parallel for
			for (i = ((44100 * 4) - 1) / 2;i < (44100 * 4);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 4));
			}
			break;
		case (48000 * 4):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 4) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 4));
			}
			#pragma omp parallel for
			for (i = ((48000 * 4) - 1) / 2;i < (48000 * 4);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 4));
			}
			break;
		case (44100 * 8):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 8) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 8));
			}
			#pragma omp parallel for
			for (i = ((44100 * 8) - 1) / 2;i < (44100 * 8);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 8));
			}
			break;
		case (48000 * 8):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 8) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 8));
			}
			#pragma omp parallel for
			for (i = ((48000 * 8) - 1) / 2;i < (48000 * 8);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 8));
			}
			break;
		case (44100 * 16):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 16) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 16));
			}
			#pragma omp parallel for
			for (i = ((44100 * 16) - 1) / 2;i < (44100 * 16);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 16));
			}
			break;
		case (48000 * 16):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 16) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 16));
			}
			#pragma omp parallel for
			for (i = ((48000 * 16) - 1) / 2;i < (48000 * 16);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 16));
			}
			break;
		case (44100 * 32):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 32) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 32));
			}
			#pragma omp parallel for
			for (i = ((44100 * 32) - 1) / 2;i < (44100 * 32);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 32));
			}
			break;
		case (48000 * 32):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 32) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 32));
			}
			#pragma omp parallel for
			for (i = ((48000 * 32) - 1) / 2;i < (48000 * 32);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 32));
			}
			break;
		case (44100 * 64):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 64) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 64));
			}
			#pragma omp parallel for
			for (i = ((44100 * 64) - 1) / 2;i < (44100 * 64);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 64));
			}
			break;
		case (48000 * 64):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 64) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 64));
			}
			#pragma omp parallel for
			for (i = ((48000 * 64) - 1) / 2;i < (48000 * 64);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 64));
			}
			break;
		case (44100 * 128):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 128) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 128));
			}
			#pragma omp parallel for
			for (i = ((44100 * 128) - 1) / 2;i < (44100 * 128);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 128));
			}
			break;
		case (48000 * 128):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 128) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 128));
			}
			#pragma omp parallel for
			for (i = ((48000 * 128) - 1) / 2;i < (48000 * 128);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 128));
			}
			break;
		case (44100 * 256):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 256) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 256));
			}
			#pragma omp parallel for
			for (i = ((44100 * 256) - 1) / 2;i < (44100 * 256);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 256));
			}
			break;
		case (48000 * 256):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 256) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 256));
			}
			#pragma omp parallel for
			for (i = ((48000 * 256) - 1) / 2;i < (48000 * 256);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 256));
			}
			break;
		default:
			for (i = 0;i < ((size) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(size));
			}
			for (i = ((size) - 1) / 2;i < (size);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(size));
			}
			break;
	}
}
//---------------------------------------------------------------------------
// Function   : cutFFTW
// Description: FFTW用カットオフ関数
// ---
//	
//
void cutFFTW(fftw_complex *fftw,long index,long size)
{
	long i;

	// 64 個ずつ
	for (i = index;i + 64 < size;i+= 64) {
		fftw[i + 0][0] = 0;
		fftw[i + 0][1] = 0;
		fftw[i + 1][0] = 0;
		fftw[i + 1][1] = 0;
		fftw[i + 2][0] = 0;
		fftw[i + 2][1] = 0;
		fftw[i + 3][0] = 0;
		fftw[i + 3][1] = 0;
		fftw[i + 4][0] = 0;
		fftw[i + 4][1] = 0;
		fftw[i + 5][0] = 0;
		fftw[i + 5][1] = 0;
		fftw[i + 6][0] = 0;
		fftw[i + 6][1] = 0;
		fftw[i + 7][0] = 0;
		fftw[i + 7][1] = 0;
		fftw[i + 8][0] = 0;
		fftw[i + 8][1] = 0;
		fftw[i + 9][0] = 0;
		fftw[i + 9][1] = 0;
		fftw[i + 10][0] = 0;
		fftw[i + 10][1] = 0;
		fftw[i + 11][0] = 0;
		fftw[i + 11][1] = 0;
		fftw[i + 12][0] = 0;
		fftw[i + 12][1] = 0;
		fftw[i + 13][0] = 0;
		fftw[i + 13][1] = 0;
		fftw[i + 14][0] = 0;
		fftw[i + 14][1] = 0;
		fftw[i + 15][0] = 0;
		fftw[i + 15][1] = 0;
		fftw[i + 16][0] = 0;
		fftw[i + 16][1] = 0;
		fftw[i + 17][0] = 0;
		fftw[i + 17][1] = 0;
		fftw[i + 18][0] = 0;
		fftw[i + 18][1] = 0;
		fftw[i + 19][0] = 0;
		fftw[i + 19][1] = 0;
		fftw[i + 20][0] = 0;
		fftw[i + 20][1] = 0;
		fftw[i + 21][0] = 0;
		fftw[i + 21][1] = 0;
		fftw[i + 22][0] = 0;
		fftw[i + 22][1] = 0;
		fftw[i + 23][0] = 0;
		fftw[i + 23][1] = 0;
		fftw[i + 24][0] = 0;
		fftw[i + 24][1] = 0;
		fftw[i + 25][0] = 0;
		fftw[i + 25][1] = 0;
		fftw[i + 26][0] = 0;
		fftw[i + 26][1] = 0;
		fftw[i + 27][0] = 0;
		fftw[i + 27][1] = 0;
		fftw[i + 28][0] = 0;
		fftw[i + 28][1] = 0;
		fftw[i + 29][0] = 0;
		fftw[i + 29][1] = 0;
		fftw[i + 30][0] = 0;
		fftw[i + 30][1] = 0;
		fftw[i + 31][0] = 0;
		fftw[i + 31][1] = 0;
		fftw[i + 32][0] = 0;
		fftw[i + 32][1] = 0;
		fftw[i + 33][0] = 0;
		fftw[i + 33][1] = 0;
		fftw[i + 34][0] = 0;
		fftw[i + 34][1] = 0;
		fftw[i + 35][0] = 0;
		fftw[i + 35][1] = 0;
		fftw[i + 36][0] = 0;
		fftw[i + 36][1] = 0;
		fftw[i + 37][0] = 0;
		fftw[i + 37][1] = 0;
		fftw[i + 38][0] = 0;
		fftw[i + 38][1] = 0;
		fftw[i + 39][0] = 0;
		fftw[i + 39][1] = 0;
		fftw[i + 40][0] = 0;
		fftw[i + 40][1] = 0;
		fftw[i + 41][0] = 0;
		fftw[i + 41][1] = 0;
		fftw[i + 42][0] = 0;
		fftw[i + 42][1] = 0;
		fftw[i + 43][0] = 0;
		fftw[i + 43][1] = 0;
		fftw[i + 44][0] = 0;
		fftw[i + 44][1] = 0;
		fftw[i + 45][0] = 0;
		fftw[i + 45][1] = 0;
		fftw[i + 46][0] = 0;
		fftw[i + 46][1] = 0;
		fftw[i + 47][0] = 0;
		fftw[i + 47][1] = 0;
		fftw[i + 48][0] = 0;
		fftw[i + 48][1] = 0;
		fftw[i + 49][0] = 0;
		fftw[i + 49][1] = 0;
		fftw[i + 50][0] = 0;
		fftw[i + 50][1] = 0;
		fftw[i + 51][0] = 0;
		fftw[i + 51][1] = 0;
		fftw[i + 52][0] = 0;
		fftw[i + 52][1] = 0;
		fftw[i + 53][0] = 0;
		fftw[i + 53][1] = 0;
		fftw[i + 54][0] = 0;
		fftw[i + 54][1] = 0;
		fftw[i + 55][0] = 0;
		fftw[i + 55][1] = 0;
		fftw[i + 56][0] = 0;
		fftw[i + 56][1] = 0;
		fftw[i + 57][0] = 0;
		fftw[i + 57][1] = 0;
		fftw[i + 58][0] = 0;
		fftw[i + 58][1] = 0;
		fftw[i + 59][0] = 0;
		fftw[i + 59][1] = 0;
		fftw[i + 60][0] = 0;
		fftw[i + 60][1] = 0;
		fftw[i + 61][0] = 0;
		fftw[i + 61][1] = 0;
		fftw[i + 62][0] = 0;
		fftw[i + 62][1] = 0;
		fftw[i + 63][0] = 0;
		fftw[i + 63][1] = 0;
	}
	// 残り
	for (;i < size;i++) {
		fftw[i + 0][0] = 0;
		fftw[i + 0][1] = 0;
	}
}
//---------------------------------------------------------------------------
// Function   : al_malloc
// Description: 16バイト境界対応malloc関数
// ---
// 返すポインタの16バイト前にmallocで確保したメモリ領域のアドレスを入れる
//
void *al_malloc(long size)
{
	void *ptr;
	void *ret;
	int align;

	ptr = malloc(size + 32);
	if (ptr) {
		ret = ptr;
		align = (int)ptr % 16;
		if (align != 0) {
			align = 16 - align;
			ret = (char *)ret + align;
		}
		*((SSIZE *)ret) = (SSIZE)ptr;

		ret = (char *)ret + 16;
	} else {
		ret = NULL;
	}
	return ret;
}
//---------------------------------------------------------------------------
// Function   : al_free
// Description: 16バイト境界対応free関数
// ---
// 
//
void *al_free(void *ptr)
{
	void *p;
	
	if (ptr) {
		p = (char *)ptr - 16;
		p = (void *)(*((SSIZE *)p));
		free(p);
	}
}
#endif
