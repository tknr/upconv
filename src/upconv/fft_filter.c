//---------------------------------------------------------------------------
/****************************************************************************/
/* upconv (C) 2007-2018 By 59414d41											*/
/*																			*/
/*																			*/
/****************************************************************************/

/*--- Log ------------------------------------------------------------------
 * Ver 0.99 <18/11/04> - upconvからソースコードを分離
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <windows.h>
#include "upconv.h"
#include "fftw3.h"
#include "fileio.h"
#include "fft_filter.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#if 0
#define	PRINT_LOG(s)	do {																	\
							FILE *log;															\
							log = fopen("d:\\fft_filter.log","a");									\
							if (log) {															\
								fprintf(log,"%s [%d] %s\n",__FUNCTION__,__LINE__,s);			\
								fclose(log);													\
							}																	\
						} while (0)
#else
#define	PRINT_LOG(s)	//
#endif
static int chkAbort(FFT_PARAM *param,int percent);

//---------------------------------------------------------------------------
// Function   : fftFilter
// Description: FFT によるフィルタ処理
// ---
//	inBuffer	:入力データバッファ
//	inSample	:入力データのサンプル数(ch毎)
//	fp_r		:入力ファイル用構造体
//	fp_w		:テンポラリファイル用構造体
//	param		:FFT変換パラメータ
//
void fftFilter(DWORD inSample,DWORD outSample,FIO *fp_r,FIO *fp_w,FFT_PARAM *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3,*mem4;
	SSIZE level,level2;
	double nx;
	SSIZE wkMemSize;
	SSIZE inSampleR,outSampleR;
	SSIZE fftSizeIn,fftSizeOut,i,j,n;
	SSIZE wr;
	double percent,per;
	SSIZE *pIn[3],*pOut[3];
	SSIZE startInSample,nSample;
	DWORD outRemain;
	fftw_complex *fftw_in[3],*fftw_out[3];
	fftw_plan fftw_p[3],fftw_ip[3];
	char s[50];

	fftw_in[0] = NULL;
	fftw_in[1] = NULL;
	fftw_in[2] = NULL;
	fftw_out[0] = NULL;
	fftw_out[1] = NULL;
	fftw_out[2] = NULL;

	fio_rewind(fp_r);

	inSampleR = param->inSampleR;
	outSampleR = param->outSampleR;

	if (param->dsd_fmt == -1) {
		fftSizeIn = inSampleR * 2;
		fftSizeOut = outSampleR * 2;
	} else {
		fftSizeIn = inSampleR / 4;
		fftSizeOut = outSampleR / 4;
	}

	wkMemSize = fftSizeIn;
	if (wkMemSize < fftSizeOut) {
		wkMemSize = fftSizeOut;
	}
	wkMemSize *= 2;

	// 入力用
	mem1 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem1 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 出力用
	mem2 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem2 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// ワーク用(別スレッド用)
	mem3 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem3 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	mem4 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem4 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 1
	fftw_in[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_in[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_out[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_out[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p[0] = fftw_plan_dft_1d(fftSizeIn,fftw_in[0],fftw_out[0],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[0] = fftw_plan_dft_1d(fftSizeOut,fftw_out[0],fftw_in[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	if (param->hi_sample_mode == 0) {
		// 2
		fftw_in[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
		if (fftw_in[1] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		fftw_out[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
		if (fftw_out[1] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		fftw_p[1] = fftw_plan_dft_1d(fftSizeIn,fftw_in[1],fftw_out[1],FFTW_FORWARD,FFTW_ESTIMATE);
		if (fftw_p[1] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		fftw_ip[1] = fftw_plan_dft_1d(fftSizeOut,fftw_out[1],fftw_in[1],FFTW_BACKWARD,FFTW_ESTIMATE);
		if (fftw_ip[1] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		// 3
		fftw_in[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
		if (fftw_in[2] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		fftw_out[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
		if (fftw_out[2] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		fftw_p[2] = fftw_plan_dft_1d(fftSizeIn,fftw_in[2],fftw_out[2],FFTW_FORWARD,FFTW_ESTIMATE);
		if (fftw_p[2] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		fftw_ip[2] = fftw_plan_dft_1d(fftSizeOut,fftw_out[2],fftw_in[2],FFTW_BACKWARD,FFTW_ESTIMATE);
		if (fftw_ip[2] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
	}

	outRemain = outSample;

	pIn[0]	= &mem1[((fftSizeIn / 2) * 0)];
	pOut[0] = &mem2[((fftSizeOut / 2) * 0)];
	pIn[1]	= &mem1[((fftSizeIn / 2) * 1)];
	pOut[1] = &mem3[((fftSizeOut / 2) * 1)];
	pIn[2]	= &mem1[((fftSizeIn / 2) * 2)];
	pOut[2] = &mem4[((fftSizeOut / 2) * 2)];
	
	per = -1;
	for (startInSample = ((fftSizeIn + (fftSizeIn / 2)) * -1);startInSample < inSample + (fftSizeIn + fftSizeIn / 2);startInSample += fftSizeIn) {
		if (startInSample >= 0 && startInSample < inSample) {
			percent = ((double)startInSample / inSample);
			percent *= 100;
			if (percent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)percent);
				fflush(stdout);
				if (chkAbort(param,percent) == 1) exit(0);
			}
			per = percent;
		}

		memset(mem1,0,wkMemSize * sizeof (SSIZE));

		if (startInSample >= 0 && startInSample + (fftSizeIn * 2) < inSample) {
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			nSample = fftSizeIn * 2;
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSizeIn * 2;
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
			nSample = fftSizeIn * 2;
		}

		memset(mem2,0,wkMemSize * sizeof (SSIZE));
		memset(mem3,0,wkMemSize * sizeof (SSIZE));
		memset(mem4,0,wkMemSize * sizeof (SSIZE));

		if (param->hi_sample_mode == 0) {
			#pragma omp parallel
			{
				#pragma omp sections
				{
					#pragma omp section
					{
						// 1
						fftFilterSub(pIn[0],pOut[0],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],param,0);
					}
					#pragma omp section
					{
						// 2
						fftFilterSub(pIn[1],pOut[1],fftw_in[1],fftw_out[1],fftw_p[1],fftw_ip[1],param,1);
					}
					#pragma omp section
					{
						// 3
						fftFilterSub(pIn[2],pOut[2],fftw_in[2],fftw_out[2],fftw_p[2],fftw_ip[2],param,2);
					}
				}
				#pragma omp for
				for (i = 0;i < wkMemSize;i++) {
					mem2[i] += mem3[i] + mem4[i];
				}
			}
		} else {
			// オーバーサンプリング(1536000)のときはメモリ消費が大きいのでfft用メモリは1つだけ確保し、並列処理はしない

			// 1
			fftFilterSub(pIn[0],pOut[0],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],param,0);
			// 2
			fftFilterSub(pIn[1],pOut[1],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],param,1);
			// 3
			fftFilterSub(pIn[2],pOut[2],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],param,2);
			#pragma omp parallel for
			for (i = 0;i < wkMemSize;i++) {
				mem2[i] += mem3[i] + mem4[i];
			}
		}

		if (startInSample + fftSizeIn / 2 >= 0) {
			// 音のレベルを調整する
			if (param->lvadj_flag == 1) {
				nx = 0.25;
				for (i = fftSizeOut / 2,n = 0;n < fftSizeOut;i++,n++) {
					mem2[i] = (SSIZE)((double)mem2[i] * nx);
				}
			}

			if (outRemain >= fftSizeOut) {
				wr = fio_write(mem2 + (fftSizeOut / 2),sizeof (SSIZE),fftSizeOut,fp_w);
				if (wr != fftSizeOut) {
					char s[100];
					sprintf(s,"%ld:fio_write(%ld,%ld)",wr,sizeof (SSIZE),fftSizeOut);
					param->err = STATUS_FILE_WRITE_ERR;
					PRINT_LOG(s);
					return;
				}
			} else {
				wr = fio_write(mem2 + (fftSizeOut / 2),sizeof (SSIZE),outRemain,fp_w);
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
	}

	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	al_free(mem4);
	
	// 1
	fftw_destroy_plan(fftw_p[0]);
	fftw_destroy_plan(fftw_ip[0]);
	fftw_free(fftw_in[0]);
	fftw_free(fftw_out[0]);

	if (param->hi_sample_mode == 0) {
		// 2
		fftw_destroy_plan(fftw_p[1]);
		fftw_destroy_plan(fftw_ip[1]);
		fftw_free(fftw_in[1]);
		fftw_free(fftw_out[1]);

		// 3
		fftw_destroy_plan(fftw_p[2]);
		fftw_destroy_plan(fftw_ip[2]);
		fftw_free(fftw_in[2]);
		fftw_free(fftw_out[2]);
	}
}
//---------------------------------------------------------------------------
// Function   : fftFilterSub
// Description: FFT によるフィルタ処理(サブ関数)
// ---
//	param		:変換パラメータ
//
void fftFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,FFT_PARAM *param,int id)
{
	SSIZE inSampleR,outSampleR;
	SSIZE wkSampleR;
	SSIZE fftSizeIn,fftSizeOut,i,ii,n;
	SSIZE cutOff;
	SSIZE hfc,lfc;
	
	double p;
	SSIZE validIndex;

	inSampleR = param->inSampleR;
	outSampleR = param->outSampleR;

	if (param->dsd_fmt == -1) {
		fftSizeIn = inSampleR * 2;
		fftSizeOut = outSampleR * 2;
	} else {
		fftSizeIn = inSampleR / 4;
		fftSizeOut = outSampleR / 4;
	}

	validIndex = ((double)(fftSizeOut * 2) / param->outSampleR) * (inSampleR / 2);
if (1) {
	char s[512];
	sprintf(s,"validIndex:%ld",validIndex);
	PRINT_LOG(s);
}
	// FFT 初期設定
	PRINT_LOG("");
	copyToFFTW(fftw_in,pIn,fftSizeIn);
	PRINT_LOG("");

	windowFFTW(fftw_in,fftSizeIn);
	PRINT_LOG("");

	for (i = 0;i < fftSizeOut;i++) {
		fftw_out[i][0] = 0;
		fftw_out[i][1] = 0;
	}
	PRINT_LOG("");
	
	fftw_execute(fftw_p);

	PRINT_LOG("");

	if (param->analyze_mode == 0) {
		// eq or 高域ディザ調整
		if (param->eq_flag || param->cut_high_dither) {
			int ii;
			for (i = 1;i < validIndex && i < fftSizeOut;i++) {
				ii = (((double)i / fftSizeOut) * outSampleR);
				if (ii < 192000 && param->eq[ii] != 0) {
					fftw_out[i][0] *= param->eq[ii];
					fftw_out[i][1] *= param->eq[ii];
				}
			}
		}
		// lfa (低域調整)
		if (param->lfa_flag == 1) {
			for (i = 1;i < validIndex && i < fftSizeOut;i++) {
				ii = (((double)i / fftSizeOut) * outSampleR);
				if (ii < 192000 && param->lfa_eq[ii] != 0) {
					fftw_out[i][0] *= param->lfa_eq[ii];
					fftw_out[i][1] *= param->lfa_eq[ii];
				}
			}
		}
		if (param->disable_eq == 0) {
			if (param->deEmphasis != 0) {
				adjPinkFilter(3,fftSizeOut,fftw_out,param);
			}
			if (param->lpf_flag != -1) {
				adjPinkFilter(4,fftSizeOut,fftw_out,param);
			}
		}
#if 0
		// スピーカーごとの調整
		if (param->ana == 1 && param->sp_ana == 0) {
			for (i = 1;i < validIndex;i++) {
				h = ((double)param->outSampleR / fftSizeOut) * i;
				if (param->sp_eq[h] != 0) {
					fftw_out[i][0] *= (param->sp_eq[h] * 0.50);
					fftw_out[i][1] *= (param->sp_eq[h] * 0.50);
				}
			}
		}
#endif
	}

	// 高域削除
	if (inSampleR <= outSampleR) {
		wkSampleR = inSampleR;
	} else {
		wkSampleR = outSampleR;
	}
	hfc = wkSampleR / 2;
	if (param->analyze_mode == 0 && param->hfc != -1) {
		hfc = param->hfc;
	}
	if (hfc > wkSampleR / 2) {
		hfc = wkSampleR / 2;
	}
	lfc = -1;
	if (param->analyze_mode == 0 && param->lfc != -1) {
		if (param->lfc + 1000 < hfc) {
			lfc = param->lfc;
		}
	}

	// 周波数補正用データ作成
	if (param->analyze_mode == 1 && id == 1) {
		for (i = 1;i < validIndex;i++) {
			p = (fftw_out[i][0] * fftw_out[i][0]) + (fftw_out[i][1] * fftw_out[i][1]);
			if (p > 0) {
				p = sqrt(p);
			}

			ii = (((double)i / fftSizeOut) * outSampleR);
			if (ii <= 19200) {
				if (param->eq_ref_max[ii] < p) {
					param->eq_ref_max[ii] = p;
				}
				param->eq_ref_avg[ii] += p;
			}
		}
		param->eq_ref_count++;
	}

	// カットオフ(hfc)
	cutOff = ((double)fftSizeOut / outSampleR) * hfc;
	cutFFTW(fftw_out,cutOff,fftSizeOut);

	if (param->analyze_mode == 0 && lfc != -1) {
		// カットオフ(lfc)
		cutOff = ((double)fftSizeOut / outSampleR) * lfc;
		for (i = 1;i < cutOff;i++) {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
	}

	// 半分のデータを復元
	PRINT_LOG("");

	for (i = 1;i < fftSizeOut / 2;i++) {
		fftw_out[fftSizeOut - i][0] = fftw_out[i][0];
		fftw_out[fftSizeOut - i][1] = fftw_out[i][1] * -1;
	}

	// invert FFT
	fftw_execute(fftw_ip);
	PRINT_LOG("");

	// 出力
	for (i = 0;i < fftSizeOut;i++) {
		pOut[i] = (SSIZE)(fftw_in[i][0] / fftSizeOut);
	}
	PRINT_LOG("");
}

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
	long i;

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
		case (4096 * 32):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 32) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 32));
			}
			#pragma omp parallel for
			for (i = ((4096 * 32) - 1) / 2;i < (4096 * 32);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 32));
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
void *al_malloc(SSIZE size)
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
//---------------------------------------------------------------------------
static int chkAbort(FFT_PARAM *param,int percent)
{
	FILE *fp;
	int diff = param->abort_percent - percent;
	
	param->abort_percent = percent;
	if (diff == 0 || diff > 10 || diff < -10 || percent == 0 || percent == 100) {
		if (param->abort_filename != NULL) {
			fp = fopen(param->abort_filename,"r");
			if (fp) {
				fclose(fp);
				return 1;
			}
		}
	}
	return 0;
}


