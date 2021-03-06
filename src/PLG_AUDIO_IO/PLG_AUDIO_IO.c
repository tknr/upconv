// ---------------------------------------------------------------------------
/** ************************************************************************* */
/* PLG_AUDIO_IO (C) 2008-2012 By 59414d41 */
/* */
/* */
/** ************************************************************************* */

/* --- Log -------------------------------------------------------------------
 * Ver 0.10 <09/06/02> - VC6用のDLLとして作成
 * Ver 0.50 <10/10/23> - MP3,FLAC,WAVPACKファイル対応(情報取得のみ対応)
 *					   - BWF 対応(保存のみ)
 *					   - RF64 対応(保存のみ)
 * Ver 0.60 <10/11/25> - Multi Channel 対応
 * Ver 0.61 <10/12/30> - Extensible 対応
 * Ver 0.70 <11/07/24> - APIを整理
 * Ver 0.80 <12/02/11> - DLL をやめて、各プログラムへリンクするように修正
 *						 mmsystemを使用しないように修正
 *						 RF64 ファイルの読み込みに対応
 * Ver 1.20 <19/11/08> - ID3タグの保存に対応(FLAC)
 */

// ---------------------------------------------------------------------------

#define _LARGEFILE_SOURCE				// for GCC
#define _FILE_OFFSET_BITS 64			// for GCC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <time.h>
#include "./../upconv/upconv.h"
#include "./PLG_AudioIO.h"
#ifndef _Windows
#include <iconv.h>
#include <sys/fcntl.h>
#include <sys/stat.h>
#else
#include <fcntl.h>
#include <sys\stat.h>
#include <io.h>
#include<sys\types.h>
#endif
#if 0
#define	PRINT_LOG(s)	do {																	\
							FILE *log;															\
							log = fopen("d:\\plg_audio.log","a");								\
							if (log) {															\
								fprintf(log,"%s [%d] %s\n","PLG",__LINE__,s);					\
								fclose(log);													\
							}																	\
						} while (0)
#else
#define	PRINT_LOG(s)	//
#endif

unsigned char tempBuffer[1 * 1024 * 1024L];

#ifndef TRUE
#define TRUE	(1)
#endif
#ifndef FALSE
#define FALSE	(0)
#endif

// Wav ファイル
#pragma pack(push, 1)
typedef struct {
	char id[4];
	DWORD size;
	unsigned char data[1];
}WAVE_HEADER;

// Wave file metadata(LIST)
typedef struct {
	BYTE h_list[4];
	DWORD s_list;
	BYTE h_info[4];
	BYTE h_isrc[4];
	DWORD s_isrc;
	BYTE *p_isrc;
	BYTE h_icrd[4];
	DWORD s_icrd;
	BYTE *p_icrd;
	BYTE h_iprd[4];
	DWORD s_iprd;
	BYTE *p_iprd;
	BYTE h_ieng[4];
	DWORD s_ieng;
	BYTE *p_ieng;
	BYTE h_inam[4];
	DWORD s_inam;
	BYTE *p_inam;
	BYTE h_icop[4];
	DWORD s_icop;
	BYTE *p_icop;
	BYTE h_ignr[4];
	DWORD s_ignr;
	BYTE *p_ignr;
	BYTE h_isft[4];
	DWORD s_isft;
	BYTE *p_isft;
	BYTE h_icmt[4];
	DWORD s_icmt;
	BYTE *p_icmt;
	BYTE h_iart[4];
	DWORD s_iart;
	BYTE *p_iart;
}META_LIST;

// Flac Picture
typedef struct {
	DWORD	type;
	DWORD	mime_size;
	BYTE	*mime_data;
	DWORD	desc_size;
	BYTE	*desc_data;
	DWORD	width;
	DWORD	height;
	DWORD	depth;
	DWORD	color;
	DWORD	pic_data_size;
	BYTE	*pic_data;
} FLAC_PIC;

// ID3 Tag,LIST(MP3,Wavに入っているタグはそのまま出力)
// Flac,DSDはWavのListとID3タグ情報を出力する(Shift-Jisは対応しない)
typedef struct {
	BYTE	id3_tag[3];
	BYTE	id3_ver[2];
	BYTE	id3_flag;
	BYTE	id3_size[4];
} ID3_H;
typedef struct {
	BYTE	talb_id[4];
	DWORD   talb_size;
	BYTE	*talb_data;
	BYTE	tpe1_id[4];
	DWORD   tpe1_size;
	BYTE	*tpe1_data;
	BYTE	tpe2_id[4];
	DWORD	tpe2_size;
	BYTE	*tpe2_data;
	BYTE	comm_id[4];
	DWORD	comm_size;
	BYTE	*comm_data;
	BYTE	tpos_id[4];
	DWORD	tpos_size;
	BYTE	*tpos_data;
	BYTE	tcon_id[4];
	DWORD	tcon_size;
	BYTE	*tcon_data;
	BYTE	tit2_id[4];
	DWORD	tit2_size;
	BYTE	*tit2_data;
	BYTE	trck_id[4];
	DWORD	trck_size;
	BYTE	*trck_data;
	BYTE	tyer_id[4];
	DWORD	tyer_size;
	BYTE	*tyer_data;
	BYTE	apic_id[4];
	DWORD	apic_size;
	BYTE	*apic_data;
} ID3_D;
#pragma pack(pop)

// サンプルを処理するデータ型
#define SSIZE	signed __int64
#define UI64	unsigned __int64

//
// DSF ファイルフォーマット仕様書を参照
#pragma pack(push, 1)
typedef struct {
	char id[4];
	UI64 chunk_size;
	UI64 file_size;
	UI64 ptr_meta;
}DSF;

typedef struct {
	char id[4];
	UI64 chunk_size;
	DWORD fmt_version;
	DWORD fmt_id;
	DWORD channel_type;
	DWORD channel_count;
	DWORD sampling;
	DWORD sample_bit_count;
	UI64 sample_count;
	DWORD block_size;
	DWORD reserved;
}DSF_FMT;

typedef struct {
	char id[4];
	UI64 chunk_size;
}DSF_DATA;
#pragma pack(pop)

int infoWavFile(char*filename, SOUNDFMT * pFmt, DWORD * pInSample,FILEINFO * pFileinfo);
int infoRF64File(char*filename, SOUNDFMT * pFmt, DWORD * pInSample,FILEINFO * pFileinfo);
int infoMP3File(char*filename, SOUNDFMT * pFmt, DWORD * pInSample,FILEINFO * pFileinfo);
int infoFlacFile(char*filename, SOUNDFMT * pFmt, DWORD * pInSample,FILEINFO * pFileinfo);
int infoWavPackFile(char*filename, SOUNDFMT * pFmt, DWORD * pInSample,FILEINFO * pFileinfo);
int infoDsfFile(char*filename, SOUNDFMT * pFmt, DWORD * pInSample,FILEINFO * pFileinfo);
int findWavChunk(FILE * fp, char*id, long*size);
int findWavChunk_N(FILE * fp,int n,char *tag,long*size);
static __int64 ftell_local(FILE * fp);
static int fseek_local(FILE * fp, __int64 offset, int origin);
char *utf8_to_utf16(char *u8str,int *u16_size);
char *utf8_to_sjis(char *u8str,int *sjis_size);
static int  u16_strlen(char *u8_str);
int infoWavFile(char *filename, SOUNDFMT *pFmt, DWORD *pInSample,FILEINFO *pFileinfo)
/*
 * WavFile の情報取得(4G以下のwavファイルのみ対応)
 */ 
{
	WAVE_HEADER *wh;
	WFEX wfex;
	__int64 seekPtr;
	DWORD pcmSize;
	DWORD chunkSize;
	char subTypePCM[16] = {
		0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xaa,
		0x00, 0x38, 0x9b, 0x71
	};
	char subTypeFloat[16] = {
		0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xaa,
		0x00, 0x38, 0x9b, 0x71
	};
	FILE *fp;
	long rs;
	int fmt_chunk, data_chunk;

	int fileIsWAV = FALSE;
	fmt_chunk = 0;
	data_chunk = 0;
	wh = (WAVE_HEADER*)tempBuffer;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// WAV
		rs = fread(tempBuffer, 1, 12, fp);
		if (rs != 12) {
			break;
		}
		if (!(memcmp(tempBuffer, "RIFF", 4) == 0 || memcmp(tempBuffer, "riff",4) == 0)) {
			break;
		}
		if (!(memcmp(tempBuffer + 8, "WAVE", 4) == 0 || memcmp(tempBuffer + 8,"wave", 4) == 0)) {
			break;
		}
//		if (wh->size > 2147483647 || wh->size + 8 > 2147483647) {
		if (wh->size > 4200000000) {
			// size over
			break;
		}

		// fmt チャンクのサーチ
		if (findWavChunk(fp, "fmt ", &chunkSize)) {
			fmt_chunk = 1;
			rs = fread(&wfex, 1, sizeof(WFEX), fp);
			if (rs < 16) {
				break;
			}
			if (wfex.Format.wFormatTag == WF_PCM) {
				if ((wfex.Format.wBitsPerSample == 16 || wfex.Format.wBitsPerSample == 24) && wfex.Format.nChannels <= 2) {
					pFmt->sample = wfex.Format.nSamplesPerSec;
					pFmt->channel = (unsigned char)wfex.Format.nChannels;
					pFmt->bitsPerSample = (unsigned char)wfex.Format.wBitsPerSample;
					fileIsWAV = TRUE;
				}
			} else if (wfex.Format.wFormatTag == WF_IEEE_FLOAT) {
				if ((wfex.Format.wBitsPerSample == 32 || wfex.Format.wBitsPerSample == 64) && wfex.Format.nChannels <= 2) {
					pFmt->sample = wfex.Format.nSamplesPerSec;
					pFmt->channel = (unsigned char)wfex.Format.nChannels;
					pFmt->bitsPerSample = (unsigned char)wfex.Format.wBitsPerSample;
					fileIsWAV = TRUE;
				}
			} else if (wfex.Format.wFormatTag == 0xFFFE && chunkSize >= sizeof(WFEX)) {
				if ((memcmp(wfex.subFormat, subTypePCM,16) == 0 && wfex.Format.wBitsPerSample == 16 || wfex.Format.wBitsPerSample == 24) ||
					(memcmp(wfex.subFormat, subTypeFloat,16) == 0 && wfex.Format.wBitsPerSample == 32 || wfex.Format.wBitsPerSample == 64)) {
					if (wfex.Format.nChannels <= 6) {
						pFmt->sample = wfex.Format.nSamplesPerSec;
						pFmt->channel = (unsigned char)wfex.Format.nChannels;
						pFmt->bitsPerSample = (unsigned char)wfex.Format.wBitsPerSample;
						fileIsWAV = TRUE;
					}
				}
			}
		}
		if (findWavChunk(fp, "data", &chunkSize)) {
			data_chunk = 1;
			seekPtr = ftell_local(fp);
			if (pFileinfo != NULL) {
				pFileinfo->dataOffset = seekPtr;
				pFileinfo->dataSize = wh->size;
			}
			pcmSize = wh->size;
			*pInSample = (pcmSize / wfex.Format.nChannels) / (wfex.Format.wBitsPerSample / 8);
		}
		if (fmt_chunk && data_chunk) {
			pFmt->fmt[0] = 'w';
			pFmt->fmt[1] = 'a';
			pFmt->fmt[2] = 'v';
			pFmt->fmt[3] = '\0';
			fileIsWAV = TRUE;
			fclose(fp);
			fp = NULL;
		} else {
			fileIsWAV = FALSE;
			fclose(fp);
			fp = NULL;
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsWAV;
}

int infoRF64File(char *filename, SOUNDFMT *pFmt, DWORD *pInSample,FILEINFO *pFileinfo)
/*
 * RF64ファイルの情報取得
 */
{
	WAVE_HEADER *wh;
	RF64_DS64 ds64;
	RF64_FMT ds64_fmt;
	__int64 seekPtr;
	__int64 pcmSize;
	DWORD chunkSize;
	FILE *fp;
	long rs;
	int fmt_chunk, data_chunk, ds64_chunk;

	int fileIsRF64 = FALSE;
	fmt_chunk = 0;
	ds64_chunk = 0;
	data_chunk = 0;
	pcmSize = 0;
	wh = (WAVE_HEADER*)tempBuffer;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// RF64
		rs = fread(tempBuffer, 1, 12, fp);
		if (rs != 12) {
			break;
		}
		if (!(memcmp(tempBuffer, "RF64", 4) == 0 || memcmp(tempBuffer, "rf64",4) == 0)) {
			break;
		}
		if (!(memcmp(tempBuffer + 8, "WAVE", 4) == 0 || memcmp(tempBuffer + 8,"wave", 4) == 0)) {
			break;
		}
		if (wh->size != 0xFFFFFFFF) {
			break;
		}
		// ds64 チャンクのサーチ
		if (findWavChunk(fp, "ds64", &chunkSize)) {
			if (fseek_local(fp, -8, SEEK_CUR)) {
				break;
			}
			ds64_chunk = 1;
			rs = fread(&ds64, sizeof(RF64_DS64), 1, fp);
			if (rs != 1) {
				break;
			}
			if (ds64.riffSizeLow == 0 && ds64.riffSizeHigh == 0) {
				break;
			}
			if (ds64.dataSizeLow == 0 && ds64.dataSizeHigh == 0) {
				break;
			}
			if (ds64.sampleCountLow == 0 && ds64.sampleCountHigh == 0) {
				break;
			}
			pcmSize = ((__int64)(ds64.dataSizeHigh) << 32) + (__int64)(ds64.dataSizeLow);
			if (pFileinfo != NULL) {
				pFileinfo->dataSize = pcmSize;
			}
		}
		// fmt チャンクのサーチ
		if (findWavChunk(fp, "fmt ", &chunkSize)) {
			if (fseek_local(fp, -8, SEEK_CUR)) {
				break;
			}
			fmt_chunk = 1;
			rs = fread(&ds64_fmt, sizeof(RF64_FMT), 1, fp);
			if (rs != 1) {
				break;
			}
			if (ds64_fmt.formatType == WF_PCM) {
				if((ds64_fmt.bitsPerSample == 16 || ds64_fmt.bitsPerSample == 24) && ds64_fmt.channelCount <= 2) {
					pFmt->sample = ds64_fmt.sampleRate;
					pFmt->channel = (unsigned char)ds64_fmt.channelCount;
					pFmt->bitsPerSample = (unsigned char)ds64_fmt.bitsPerSample;
					fileIsRF64 = TRUE;
				}
			}
			else if (ds64_fmt.formatType == WF_IEEE_FLOAT) {
				if((ds64_fmt.bitsPerSample == 32 || ds64_fmt.bitsPerSample == 64) && ds64_fmt.channelCount <= 2) {
					pFmt->sample = ds64_fmt.sampleRate;
					pFmt->channel = (unsigned char)ds64_fmt.channelCount;
					pFmt->bitsPerSample = (unsigned char)ds64_fmt.bitsPerSample;
					fileIsRF64 = TRUE;
				}
			}
		}
		if (findWavChunk(fp, "data", &chunkSize)) {
			data_chunk = 1;
			seekPtr = ftell_local(fp);
			if (pFileinfo != NULL) {
				pFileinfo->dataOffset = seekPtr;
			}
		}
		if (fmt_chunk && data_chunk && ds64_chunk) {
			pFmt->fmt[0] = 'r';
			pFmt->fmt[1] = 'f';
			pFmt->fmt[2] = '6';
			pFmt->fmt[3] = '4';
			pFmt->fmt[4] = '\0';
			PRINT_LOG("");
			*pInSample = ((pcmSize / pFmt->channel) / (pFmt->bitsPerSample / 8)) <= 0xFFFFFFFF ? ((pcmSize / pFmt->channel) / (pFmt->bitsPerSample / 8)) : 0;
			if (*pInSample > 0) {
				// サンプル数を32bitに制限する(本ソフトの仕様)
				fileIsRF64 = TRUE;
			}
			fclose(fp);
			fp = NULL;
		} else {
			fileIsRF64 = FALSE;
			fclose(fp);
			fp = NULL;
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsRF64;
}

int infoMP3File(char *filename, SOUNDFMT *pFmt, DWORD *pInSample,
	FILEINFO *pFileinfo)
/*
 * MP3 Fileの情報取得
 */ {
	FILE *fp;
	long rs;
	__int64 seekPtr;
	int fileIsMP3 = FALSE;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// MP3
		rs = fread(tempBuffer, 1, 10, fp);
		if (rs != 10) {
			break;
		}
		if (tempBuffer[0] == 'I' && tempBuffer[1] == 'D' && tempBuffer[2] == '3') {
			// ID3 v2 Tag
			seekPtr = (tempBuffer[6] & 0x7F);
			seekPtr <<= 7;
			seekPtr |= (tempBuffer[7] & 0x7F);
			seekPtr <<= 7;
			seekPtr |= (tempBuffer[8] & 0x7F);
			seekPtr <<= 7;
			seekPtr |= (tempBuffer[9] & 0x7F);
			seekPtr += 10;
			if (fseek_local(fp, seekPtr, SEEK_SET)) {
				break;
			}
			rs = fread(tempBuffer, 1, 4, fp);
			if (rs != 4) {
				break;
			}
		}
		if (tempBuffer[0] == 0xFF && (tempBuffer[1] & 0xF0) == 0xF0 && (tempBuffer[1] & 0x08) == 0x08) {
			int layer = (tempBuffer[1] >> 1) & 0x03;
			int sf = (tempBuffer[2] >> 2) & 0x03;
			long st[3] = {
				44100, 48000, 32000
			};
			int ch = (tempBuffer[3] >> 6) & 0x03;
			if (sf >= 0 && sf < 3) {
				*pInSample = 0;
				pFmt->sample = st[sf];
				pFmt->channel = 2;
				if (ch == 3) {
					pFmt->channel = 1;
				}
				pFmt->bitsPerSample = 16;
				pFmt->fmt[0] = 'm';
				pFmt->fmt[1] = 'p';
				pFmt->fmt[2] = '3';
				pFmt->fmt[3] = '\0';
				fileIsMP3 = TRUE;
				fclose(fp);
				fp = NULL;
			}
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsMP3;
}

int infoFlacFile(char *filename, SOUNDFMT *pFmt, DWORD *pInSample,FILEINFO *pFileinfo)
/*
 * Flac Fileの情報取得
 */ {
	FILE *fp;
	long rs;
	int fileIsFLAC = FALSE;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// Flac
		rs = fread(tempBuffer, 1, 10, fp);
		if (rs != 10) {
			break;
		}
		if (tempBuffer[0] == 'f' && tempBuffer[1] == 'L' && tempBuffer[2] == 'a' && tempBuffer[3] == 'C') {
			if (fseek_local(fp, 4, SEEK_SET)) {
				break;
			}
			rs = fread(tempBuffer, 1, 38, fp);
			if (rs != 38) {
				break;
			}
			if ((tempBuffer[0] & 0x7F) == 0x00) {
				int sr;
				int ch;
				int bit;
				sr = tempBuffer[14];
				sr <<= 8;
				sr |= tempBuffer[15];
				sr <<= 4;
				sr |= (tempBuffer[16] >> 4) & 0x0F;
				ch = (tempBuffer[16] >> 1) & 0x07;
				ch++;
				bit = (tempBuffer[16] & 0x01) << 4;
				bit |= (tempBuffer[17] >> 4) & 0x0F;
				bit++;
				*pInSample = 0;
				pFmt->sample = sr;
				pFmt->channel = ch;
				pFmt->bitsPerSample = bit;
				pFmt->fmt[0] = 'f';
				pFmt->fmt[1] = 'l';
				pFmt->fmt[2] = 'a';
				pFmt->fmt[3] = 'c';
				pFmt->fmt[4] = '\0';
				fileIsFLAC = TRUE;
				fclose(fp);
				fp = NULL;
			}
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsFLAC;
}

int infoWavPackFile(char *filename, SOUNDFMT *pFmt, DWORD *pInSample,
	FILEINFO *pFileinfo)
/*
 * WavPack File の情報取得
 */ {
	FILE *fp;
	long rs;
	int fileIsWV = FALSE;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// WavPack
		rs = fread(tempBuffer, 1, 10, fp);
		if (rs != 10) {
			break;
		}
		if (tempBuffer[0] == 'w' && tempBuffer[1] == 'v' && tempBuffer[2] == 'p' && tempBuffer[3] == 'k') {
			if (fseek_local(fp, 24, SEEK_SET)) {
				break;
			}
			rs = fread(tempBuffer, 1, 4, fp);
			if (rs == 4) {
				int bit[4] = {
					8, 16, 24, 32
				};
				long st[16] = {6000, 8000, 9600, 11025, 12000, 16000, 22050, 24000, 32000, 44100, 48000, 64000, 88200, 96000, 192000, 0};
				int stIndex;
				int bitIndex;
				*pInSample = 0;
				bitIndex = (tempBuffer[0] & 0x03);
				stIndex = ((tempBuffer[2] >> 7) & 0x01) | ((unsigned short)tempBuffer[3] << 1);
				if (bitIndex != 0 && stIndex >= 8 && stIndex <= 14 && stIndex != 11) {
					pFmt->sample = st[stIndex];
					pFmt->channel = ((tempBuffer[0] >> 2) & 1) == 0 ? 2 : 1;
					pFmt->bitsPerSample = bit[bitIndex];
					pFmt->fmt[0] = 'w';
					pFmt->fmt[1] = 'a';
					pFmt->fmt[2] = 'v';
					pFmt->fmt[3] = 'p';
					pFmt->fmt[4] = '\0';
					fileIsWV = TRUE;
					fclose(fp);
					fp = NULL;
				}
			}
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsWV;
}

int infoDsfFile(char *filename, SOUNDFMT *pFmt, DWORD *pInSample,
	FILEINFO *pFileinfo)
/*
 * DSF File の情報取得
 */ {
	FILE *fp;
	long rs;
	DSF dsf;
	DSF_FMT dsf_fmt;
	DSF_DATA dsf_data;
	int fileIsDSF = FALSE;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// DSF
		rs = fread(&dsf, 1, sizeof(DSF), fp);
		if (rs != sizeof(DSF)) {
			break;
		}
		if (memcmp(dsf.id, "DSD ", 4)) {
			break;
		}
		if (dsf.chunk_size < 28) {
			break;
		}
		if (dsf.file_size < (28 + 52 + 12)) {
			break;
		}

		if (fseek_local(fp, dsf.chunk_size, SEEK_SET)) {
			break;
		}
		rs = fread(&dsf_fmt, 1, sizeof(DSF_FMT), fp);
		if (rs != sizeof(DSF_FMT)) {
			break;
		}
		if (memcmp(dsf_fmt.id, "fmt ", 4)) {
			break;
		}
		if (dsf_fmt.chunk_size < 52) {
			break;
		}
		if (dsf_fmt.fmt_version != 1) {
			break;
		}
		if (dsf_fmt.fmt_id != 0) {
			break;
		}
		if (!(dsf_fmt.channel_type == 1 || dsf_fmt.channel_type == 2)) {
			break;
		}
		if (!(dsf_fmt.channel_count == 1 || dsf_fmt.channel_count == 2)) {
			break;
		}
		if (!(dsf_fmt.sampling == 2822400 || dsf_fmt.sampling == 2822400 * 2 || dsf_fmt.sampling == 2822400 * 4)) {
			break;
		}
		if (dsf_fmt.sample_bit_count != 1) {
			break;
		}
		if (dsf_fmt.sample_count < 1) {
			break;
		}
		if (dsf_fmt.block_size != 4096) {
			break;
		}
		if (fseek_local(fp, dsf.chunk_size + dsf_fmt.chunk_size, SEEK_SET)) {
			break;
		}
		rs = fread(&dsf_data, 1, sizeof(DSF_DATA), fp);
		if (rs != sizeof(DSF_DATA)) {
			break;
		}
		if (dsf_data.chunk_size <= 12) {
			break;
		}
		*pInSample = dsf_fmt.sample_count / dsf_fmt.sampling;
		pFmt->sample = dsf_fmt.sampling;
		pFmt->channel = dsf_fmt.channel_count;
		pFmt->bitsPerSample = dsf_fmt.sample_bit_count;
		pFmt->fmt[0] = 'd';
		pFmt->fmt[1] = 's';
		pFmt->fmt[2] = 'f';
		pFmt->fmt[3] = '\0';
		fileIsDSF = TRUE;
		fclose(fp);
		fp = NULL;
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsDSF;

}

/* */
// ---------------------------------------------------------------------------
/*
 * Function    : PLG_InfoAudioData
 * Description : Wave ファイルの情報を得る
 *				 サポートするフォーマットは16,20,24,32(float),64(float) のファイル(それ以外はサポートしない)
 *				 RF64ファイル
 *				 MP3ファイル
 *				 Flacファイル
 *				 Wavpackファイル
 *				 DSDファイル
 * ---
 *	filename   : ファイル名
 *	pFmt	   : 音声ファイル情報を格納する構造体のポインタ
 *	pInSample  : 1Ch あたりのサンプル数を格納する変数のポインタ
 *	pFileinfo  : 曲情報を格納する構造体のポインタ
 *
 */
int PLG_InfoAudioData(char *filename, SOUNDFMT *pFmt, DWORD *pInSample,
	FILEINFO *pFileinfo)
/*
 * Audio Data の情報取得
 */ {
	int fd;
	struct stat st;
	struct tm *s_time;
	int retValue;

	memset(pFmt, 0, sizeof(SOUNDFMT));
	if (pFileinfo != NULL) {
		memset(pFileinfo, 0, sizeof(FILEINFO));
	}
	fd = open(filename, O_RDONLY);
	if (fd >= 0) {
		if (!fstat(fd, &st)) {
			s_time = localtime(&st.st_mtime);
			sprintf(pFmt->date, "%4d-%02d-%02d", 1900 + s_time->tm_year,s_time->tm_mon, s_time->tm_mday);
			sprintf(pFmt->time, "%02d:%02d:%02d", s_time->tm_hour,s_time->tm_min, s_time->tm_sec);
		}
		close(fd);
	}
	retValue = STATUS_FILE_READ_ERR;

	do {
		retValue = STATUS_UNKNOWN_FORMAT;
		if (infoWavFile(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoRF64File(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoMP3File(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoFlacFile(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoWavPackFile(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoDsfFile(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}

	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_MakeHeaderWAV
 * Description : Wav ファイルのヘッダを作成する
 * ---
 *	pInFmt		: Wave ファイル情報を格納する構造体のポインタ(元)
 *	pOutFmt 	: Wave ファイル情報を格納する構造体のポインタ(先)
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	size		: bufferのサイズ
 *
 */
int PLG_MakeHeaderWAV(SOUNDFMT *pInFmt, SOUNDFMT *pOutFmt, char *buffer,
	long size, long *header_size)
/*
 * Wav ヘッダ作成
 */ {
	int retValue;
	char subTypePCM[16] = {0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xaa,0x00, 0x38, 0x9b, 0x71};
	char subTypeFloat[16] = {0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xaa,0x00, 0x38, 0x9b, 0x71};
	WAVE_HEADER *wh;
	WFE *wfe;
	WFEX *wfex;
	long ws;

	retValue = STATUS_MEM_ALLOC_ERR;

	do {
		if (pOutFmt->channel <= 2) {
			if (size <= 20+sizeof(WFE)) {
				break;
			}
		} else {
			if (size <= 20+sizeof(WFEX)) {
				break;
			}
		}
		ws = 0;
		wh = (WAVE_HEADER*)buffer;
		wh->id[0] = 'R';
		wh->id[1] = 'I';
		wh->id[2] = 'F';
		wh->id[3] = 'F';
		wh->size = 0;
		wh->data[0] = 'W';
		wh->data[1] = 'A';
		wh->data[2] = 'V';
		wh->data[3] = 'E';
		//
		wh->data[4] = 'f';
		wh->data[5] = 'm';
		wh->data[6] = 't';
		wh->data[7] = ' ';
		if (pOutFmt->channel <= 2) {
			wh->data[8] = (unsigned char)(sizeof(WFE));
			wh->data[9] = (unsigned char)(sizeof(WFE) >> 8);
			wh->data[10] = (unsigned char)(sizeof(WFE) >> 16);
			wh->data[11] = (unsigned char)(sizeof(WFE) >> 24);
			ws += 20;
			wfe = (WFE*) & buffer[ws];
			memset(wfe, 0, sizeof(WFE));
			/* Set WaveformatEx */
			if (pOutFmt->bitsPerSample == 16 || pOutFmt->bitsPerSample == 24) {
				wfe->wFormatTag = WF_PCM;
			} else if (pOutFmt->bitsPerSample == 32 || pOutFmt->bitsPerSample == 64) {
				wfe->wFormatTag = WF_IEEE_FLOAT;
			}
			wfe->nChannels = pOutFmt->channel;
			wfe->nSamplesPerSec = pOutFmt->sample;
			wfe->wBitsPerSample = pOutFmt->bitsPerSample;
			wfe->nBlockAlign = pOutFmt->channel * pOutFmt->bitsPerSample / 8;
			wfe->nAvgBytesPerSec = pOutFmt->sample * wfe->nBlockAlign;
			wfe->cbSize = 0;
			ws += sizeof(WFE);
		} else {
			wh->data[8] = (unsigned char)(sizeof(WFEX));
			wh->data[9] = (unsigned char)(sizeof(WFEX) >> 8);
			wh->data[10] = (unsigned char)(sizeof(WFEX) >> 16);
			wh->data[11] = (unsigned char)(sizeof(WFEX) >> 24);
			ws += 20;
			wfex = (WFEX*) & buffer[ws];
			memset(wfex, 0, sizeof(WFEX));
			if (pOutFmt->bitsPerSample == 16 || pOutFmt->bitsPerSample == 24) {
				memcpy(wfex->subFormat, subTypePCM, 16);
			} else if (pOutFmt->bitsPerSample == 32 || pOutFmt->bitsPerSample == 64) {
				memcpy(wfex->subFormat, subTypeFloat, 16);
			}
			wfex->Format.wFormatTag = 0xFFFE;
			wfex->Format.nChannels = pOutFmt->channel;
			wfex->Format.nSamplesPerSec = pOutFmt->sample;
			wfex->Format.wBitsPerSample = pOutFmt->bitsPerSample;
			wfex->Format.nBlockAlign = wfex->Format.nChannels * pOutFmt->bitsPerSample / 8;
			wfex->Format.nAvgBytesPerSec = pOutFmt->sample * wfex->Format.nBlockAlign;
			wfex->Format.cbSize = 22;
			wfex->Samples.wValidBitsPerSample = pOutFmt->bitsPerSample;

			ws += sizeof(WFEX);
		}
		wh = (WAVE_HEADER*) & buffer[ws];
		wh->id[0] = 'd';
		wh->id[1] = 'a';
		wh->id[2] = 't';
		wh->id[3] = 'a';
		wh->size = 0;
		ws += 8;
		*header_size = ws;
		retValue = STATUS_SUCCESS;
	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_UpdateHeaderWAV
 * Description : Wav ファイルのヘッダを更新する
 * ---
 *	filesize	: ファイルサイズ
 *	datasize	: 音声データサイズ
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	header_size	: bufferのサイズ
 *
 */
int PLG_UpdateHeaderWAV(SOUNDFMT *pOutFmt, __int64 filesize, __int64 datasize,char *buffer, long header_size)
/*
 * Wav ヘッダ更新
 */
{
	int retValue;
	WAVE_HEADER *wh;
	WFE *wfe;
	WFEX *wfex;
	long ws;

	retValue = STATUS_MEM_ALLOC_ERR;
	do {
		ws = 0;
		if (pOutFmt->channel <= 2) {
			if (header_size < 20+sizeof(WFE)) {
				break;
			}
		} else {
			if (header_size < 20+sizeof(WFEX)) {
				break;
			}
		}
		wh = (WAVE_HEADER*)buffer;
		wh->size = filesize - 8;
		//
		ws += 12 + 8;
		if (pOutFmt->channel <= 2) {
			ws += sizeof(WFE);
		} else {
			ws += sizeof(WFEX);
		}
		wh = (WAVE_HEADER*) & buffer[ws];
		wh->size = datasize;
		retValue = STATUS_SUCCESS;
	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_MakeHeaderRF64
 * Description : RF64 ファイルのヘッダを作成する
 * ---
 *	pInFmt		: Wave ファイル情報を格納する構造体のポインタ(元)
 *	pOutFmt 	: Wave ファイル情報を格納する構造体のポインタ(先)
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	size		: bufferのサイズ
 *
 */
int PLG_MakeHeaderRF64(SOUNDFMT *pInFmt, SOUNDFMT *pOutFmt, char *buffer,long size, long *header_size)
/*
 * Wav ヘッダ作成
 */
{
	int retValue;
	RF64 *rf64;
	RF64_DS64 *ds64;
	RF64_FMT *fmt;
	RF64_DATA *data;
	__int64 ws;

	retValue = STATUS_MEM_ALLOC_ERR;

	do {
		if (size < 12+sizeof(RF64_DS64)+sizeof(RF64_FMT) + (sizeof(RF64_DATA) - 1)) {
			break;
		}
		ws = 0;

		rf64 = (RF64*)buffer;
		/* RF64 */
		rf64->chunkId[0] = 'R';
		rf64->chunkId[1] = 'F';
		rf64->chunkId[2] = '6';
		rf64->chunkId[3] = '4';
		rf64->chunkSize = 0xFFFFFFFF;
		rf64->type[0] = 'W';
		rf64->type[1] = 'A';
		rf64->type[2] = 'V';
		rf64->type[3] = 'E';

		/* ds64 */
		ws += sizeof(RF64);
		ds64 = (RF64_DS64*) & buffer[ws];
		memset(ds64, 0, sizeof(RF64_DS64));
		ds64->chunkId[0] = 'd';
		ds64->chunkId[1] = 's';
		ds64->chunkId[2] = '6';
		ds64->chunkId[3] = '4';
		ds64->chunkSize = sizeof(RF64_DS64) - 8;

		/* fmt */
		ws += sizeof(RF64_DS64);
		fmt = (RF64_FMT*) & buffer[ws];
		fmt->chunkId[0] = 'f';
		fmt->chunkId[1] = 'm';
		fmt->chunkId[2] = 't';
		fmt->chunkId[3] = ' ';
		fmt->chunkSize = sizeof(RF64_FMT) - 8;
		if (pOutFmt->bitsPerSample == 16 || pOutFmt->bitsPerSample == 24) {
			fmt->formatType = WF_PCM;
		}
		else if (pOutFmt->bitsPerSample == 32 || pOutFmt->bitsPerSample == 64) {
			fmt->formatType = WF_IEEE_FLOAT;
		}
		fmt->channelCount = pOutFmt->channel;
		fmt->sampleRate = pOutFmt->sample;
		fmt->bitsPerSample = pOutFmt->bitsPerSample;
		fmt->blockAlignment = pOutFmt->channel * pOutFmt->bitsPerSample / 8;
		fmt->bytesPerSecond = pOutFmt->sample * fmt->blockAlignment;

		/* data */
		ws += sizeof(RF64_FMT);
		data = (RF64_DATA*) & buffer[ws];
		data->chunkId[0] = 'd';
		data->chunkId[1] = 'a';
		data->chunkId[2] = 't';
		data->chunkId[3] = 'a';
		data->chunkSize = 0xFFFFFFFF;
		ws += 8;
		*header_size = ws;
		retValue = STATUS_SUCCESS;
	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_UpdateHeaderRF64
 * Description : Wav ファイルのヘッダを更新する
 * ---
 *	filesize	: ファイルサイズ
 *	datasize	: 音声データサイズ
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	header_size	: bufferのサイズ
 *
 */
int PLG_UpdateHeaderRF64(SOUNDFMT *pOutFmt, __int64 filesize, __int64 datasize,char *buffer, long header_size)
/*
 * Wav ヘッダ更新
 */
{
	int retValue;
	RF64_DS64 *ds64;
	DWORD dwLow, dwHigh;
	__int64 ws;
	__int64 samplecount;

	retValue = STATUS_MEM_ALLOC_ERR;

	do {
		if (header_size < 12+sizeof(RF64_DS64)+sizeof(RF64_FMT) + (sizeof(RF64_DATA) - 1)) {
			break;
		}
		ws = 0;

		/* ds64 */
		ws += sizeof(RF64);
		ds64 = (RF64_DS64*) & buffer[ws];
		dwLow = (DWORD)(filesize - 8);
		dwHigh = (DWORD)((filesize - 8) >> 32);
		ds64->riffSizeLow = dwLow;
		ds64->riffSizeHigh = dwHigh;
		dwLow = (DWORD)datasize;
		dwHigh = (DWORD)(datasize >> 32);
		ds64->dataSizeLow = dwLow;
		ds64->dataSizeHigh = dwHigh;
		samplecount = datasize;
		samplecount /= (pOutFmt->channel * (pOutFmt->bitsPerSample / 8));
		ds64->sampleCountLow = (DWORD)(samplecount);
		ds64->sampleCountHigh = (DWORD)((samplecount >> 32));

		retValue = STATUS_SUCCESS;
	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_MakeHeaderDSD
 * Description : Dsd ファイルのヘッダを作成する
 * ---
 *	pInFmt		: DSD ファイル情報を格納する構造体のポインタ(元)
 *	pOutFmt 	: DSD ファイル情報を格納する構造体のポインタ(先)
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	size		: bufferのサイズ
 *
 */
int PLG_MakeHeaderDSD(SOUNDFMT *pInFmt, SOUNDFMT *pOutFmt, char *buffer, long size, long *header_size)
/*
 * DSF File の情報取得
 */
{
	int retValue;
	DSF *dsf;
	DSF_FMT *dsf_fmt;
	DSF_DATA *dsf_data;
	__int64 ws;

	retValue = STATUS_MEM_ALLOC_ERR;

	do {
		if (size < (sizeof (DSF) + sizeof (DSF_FMT) + sizeof (DSF_DATA))) {
			break;
		}

		retValue = STATUS_DATA_ERR;
		if (!(pOutFmt->channel == 1 || pOutFmt->channel == 2)) {
			break;
		}
		if (!(pOutFmt->sample == 2822400 || pOutFmt->sample == 2822400 * 2 || pOutFmt->sample == 2822400 * 4)) {
			break;
		}

		ws = 0;
		dsf = (DSF *)buffer;
		memcpy(dsf->id, "DSD ", 4);
		dsf->chunk_size = 28;
		dsf->file_size = sizeof (DSF) + sizeof (DSF_FMT) + sizeof (DSF_DATA);
		dsf->ptr_meta  = dsf->file_size;
		ws += sizeof (DSF);

		dsf_fmt = (DSF_FMT *)&buffer[ws];

		memcpy(dsf_fmt->id,"fmt ", 4);
		dsf_fmt->chunk_size = sizeof (DSF_FMT);
		dsf_fmt->fmt_version = 1;
		dsf_fmt->fmt_id = 0;
		dsf_fmt->channel_type  = pOutFmt->channel;
		dsf_fmt->channel_count = pOutFmt->channel;
		dsf_fmt->sampling = pOutFmt->sample;
		dsf_fmt->sample_bit_count = 1;
		dsf_fmt->sample_count = 0;
		dsf_fmt->block_size = 4096;
		ws += sizeof (DSF_FMT);

		dsf_data = (DSF_DATA *)&buffer[ws];
		memcpy(dsf_data->id,"data",4);
		dsf_data->chunk_size = sizeof (DSF_DATA);
		ws += sizeof (DSF_DATA);
		*header_size = ws;
		retValue = STATUS_SUCCESS;
	} while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_UpdateHeaderDSD
 * Description : DSD ファイルのヘッダを更新する
 * ---
 *	filesize	: ファイルサイズ
 *	datasize	: 音声データサイズ
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	header_size	: bufferのサイズ
 *
 */
int PLG_UpdateHeaderDSD(SOUNDFMT *pOutFmt, __int64 filesize, __int64 datasize,__int64 sample_count,char *buffer, long header_size)
/*
 * DSF ヘッダ更新
 */
{
	int retValue;
	DSF *dsf;
	DSF_FMT *dsf_fmt;
	DSF_DATA *dsf_data;
	long ws;
PRINT_LOG("");
	retValue = STATUS_MEM_ALLOC_ERR;
	do {
		ws = 0;
		if (pOutFmt->channel <= 2) {
			if (header_size < 20+sizeof(WFE)) {
				break;
			}
		} else {
			if (header_size < 20+sizeof(WFEX)) {
				break;
			}
		}
PRINT_LOG("");
		dsf = (DSF*)&buffer[0];
		dsf->file_size = filesize;
		dsf->ptr_meta  = filesize;
		ws += sizeof (DSF);
PRINT_LOG("");

		dsf_fmt = (DSF_FMT *)&buffer[ws];
		dsf_fmt->sample_count = sample_count;

		ws += sizeof (DSF_FMT);
PRINT_LOG("");

		dsf_data = (DSF_DATA *)&buffer[ws];
PRINT_LOG("");
		dsf_data->chunk_size = sizeof (DSF_DATA) + datasize;
PRINT_LOG("");
		retValue = STATUS_SUCCESS;
	}
	while (0);
PRINT_LOG("");

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_GetExtraChunkFor
 * Description : Wav,mp3,flac ファイルのfmt,data以外のチャンクを取得する
 * ---
 *	filesize	: ファイルサイズ
 *	nChunk		: 取得するチャンクの番号
 *	infoChunk	: 取得する情報を入れるバッファ
 *	infoSize	: 返すサイズ
 *
 */
int PLG_GetExtraChunk(char *filename, int nChunk, char **infoChunk,long *infoSize)
{
	int retValue;
	char *p1, *p2;
	BYTE header[12];
	BYTE *ovcom;
	long ovcom_count;
	FILE *fp;
	SOUNDFMT inFmt;
	DWORD inSample;
	FILEINFO fileinfo;
	WAVE_HEADER *wh;
	META_LIST meta_list;
	ID3_H id3_h;
	ID3_D id3_d;
	FLAC_PIC flac_pic;
	DWORD rs;
	DWORD sizeLow,sizeHigh;
	__int64 dataSize;
	__int64 seekPtr;
	DWORD wk_long;
	int i, next, end;
	int err;
	int flg_datachunk = 0;
	int index;
	BYTE *p_u16;
	BYTE *p_sjis;
	int  u16_size;
	int  sjis_size;
	BYTE b0,b1,b2,b3;

	retValue = STATUS_FILE_READ_ERR;

	memset(&id3_h,0,sizeof(ID3_H));
	memset(&id3_d,0,sizeof(ID3_D));
	memset(&flac_pic,0,sizeof(FLAC_PIC));

	if (*infoChunk) {
		free(*infoChunk);
		*infoChunk = NULL;
	}
	*infoSize = 0;
	dataSize  = 0;
	sizeLow = sizeHigh = 0;

	wh = (WAVE_HEADER*)header;
	ovcom = NULL;
	p_u16 = p_sjis = NULL;
	do {
		memset(&meta_list, 0, sizeof(META_LIST));

		fp = NULL;
		if (infoWavFile(filename, &inFmt, &inSample, &fileinfo) || infoRF64File(filename, &inFmt, &inSample, &fileinfo)) {
			PRINT_LOG("infoWavFile");
			fp = fopen(filename, "rb");
			if (fp == NULL) {
				break;
			}
			// ヘッダ読み込み
			rs = fread(wh, 1, 12, fp);
			if (rs != 12) {
				break;
			}
			seekPtr = 12;
			// 指定チャンクまでスキップ
			next = 0;
			for (i = 0; i < nChunk; ) {
				err = 0;
				if (next) {
					if (flg_datachunk == 0) {
						if ((wh->size & 1)) {
							seekPtr += wh->size + 9;
						} else {
							seekPtr += wh->size + 8;
						}
					} else {
						flg_datachunk = 0;
						seekPtr += dataSize + 8;
					}
					next = 0;
				}
				if (fseek_local(fp, seekPtr, SEEK_SET)) {
					err = 1;
					break;
				}
				rs = fread(wh, 1, 8, fp);
				if (rs != 8) {
					err = 1;
					break;
				}
				if (!((wh->id[0] == 'f' && wh->id[1] == 'm' && wh->id[2] == 't' && wh->id[3] == ' ') ||
					(wh->id[0] == 'd' && wh->id[1] == 'a' && wh->id[2] == 't' && wh->id[3] == 'a')   ||
					(wh->id[0] == 'd' && wh->id[1] == 's' && wh->id[2] == '6' && wh->id[3] == '4')  ||
					(wh->id[0] == 'J' && wh->id[1] == 'U' && wh->id[2] == 'N' && wh->id[3] == 'K'))) {
					i++;
				} else if (wh->id[0] == 'd' && wh->id[1] == 's' && wh->id[2] == '6' && wh->id[3] == '4') {
					rs = fread(&sizeLow,4,1,fp);
					if (rs != 1) {
						err = 1;
						break;
					}
					rs = fread(&sizeHigh,4,1,fp);
					if (rs != 1) {
						err = 1;
						break;
					}
					rs = fread(&sizeLow,4,1,fp);
					if (rs != 1) {
						err = 1;
						break;
					}
					rs = fread(&sizeHigh,4,1,fp);
					if (rs != 1) {
						err = 1;
						break;
					}
					dataSize = sizeHigh;dataSize <<= 32;
					dataSize |= sizeLow;
					if (fseek_local(fp, seekPtr, SEEK_SET)) {
						break;
					}
					err = 1;
				} else if (wh->id[0] == 'd' && wh->id[1] == 'a' && wh->id[2] == 't' && wh->id[3] == 'a') {
					flg_datachunk = 1;
					if (wh->size != 0xFFFFFFFF) {
						dataSize = wh->size;
					}
					err = 1;
				} else {
					err = 1;
				}
				next = 1;
			}
			if (err == 0 && wh->size > 0) {
				*infoChunk = (char*)malloc(wh->size + 9);
				if (*infoChunk) {
					memset(*infoChunk,0,wh->size + 9);
					if (fseek_local(fp, seekPtr, SEEK_SET) == 0) {
						rs = fread(*infoChunk, 1, wh->size + 8, fp);
						if (rs == wh->size + 8) {
							*infoSize = wh->size + 8;
						}
					}
					if (*infoSize != wh->size + 8) {
						free(*infoChunk);
						*infoSize = 0;
					}
				}
			}
		} else if (infoMP3File(filename, &inFmt, &inSample, &fileinfo) == TRUE) {
			PRINT_LOG("infoMP3File");

			if (nChunk > 1) break;

			fp = fopen(filename, "rb");
			if (fp == NULL) {
				break;
			}
			// ヘッダ読み込み
			rs = fread(header,1,10, fp);
			if (rs != 10) {
				break;
			}
			if (header[0] == 'I' && header[1] == 'D' && header[2] == '3' && header[3] == 3 && header[4] == 0) {
				b3 = header[6] & 0x7F;
				b2 = header[7] & 0x7F;
				b1 = header[8] & 0x7F;
				b0 = header[9] & 0x7F;
				wk_long  = (long)b3 << (7 * 3);
				wk_long |= (long)b2 << (7 * 2);
				wk_long |= (long)b1 << (7 * 1);
				wk_long |= b0;

				*infoSize = wk_long + 18;
				*infoChunk = malloc(*infoSize);
				if (*infoChunk != NULL) {
					p1 = *infoChunk;p1 += 8;
					memcpy(p1,header,10);p1 += 10;
					rs = fread(p1,1,wk_long,fp);
					if (rs != wk_long) {
						free(*infoChunk);
						*infoChunk = NULL;
						*infoSize = 0;
						PRINT_LOG("not wk_long");
						break;
					} else {
						p1 = *infoChunk + 18;
						if (wk_long > 0) {
							char s[100];
							// パディング領域のサイズも含むためサイズを算出する
							long id3_size,remain;
							id3_size = 0;
							remain = wk_long;
							p1 = *infoChunk + 18;
							while (remain > 0) {
								wk_long = ((long)(BYTE)p1[4] << 24);
								wk_long |= (long)(BYTE)p1[5] << 16;
								wk_long |= (long)(BYTE)p1[6] <<  8;
								wk_long |= (long)(BYTE)p1[7];
								sprintf(s,"%c%c%c%c:%08lx",p1[0],p1[1],p1[2],p1[3],wk_long);
								PRINT_LOG(s);
								if (p1[0] == 0 && p1[1] == 0 && p1[2] == 0 && p1[3] == 0 && wk_long == 0) break;
								p1 += (10 + wk_long);
								id3_size += (10 + wk_long);
								remain   -= (10 + wk_long);
							}
							wk_long = id3_size;
						}
						p1 = *infoChunk;
						p1[0] = 'i';
						p1[1] = 'd';
						p1[2] = '3';
						p1[3] = ' ';
						p1[4] = (BYTE)((wk_long + 10) >> 0);
						p1[5] = (BYTE)((wk_long + 10) >> 8);
						p1[6] = (BYTE)((wk_long + 10) >> 16);
						p1[7] = (BYTE)((wk_long + 10) >> 24);
						p1   += 8;
						p1[6] = (BYTE)(wk_long >> (7 * 3));
						p1[7] = (BYTE)(wk_long >> (7 * 2));
						p1[8] = (BYTE)(wk_long >> (7 * 1));
						p1[9] = (BYTE)((wk_long & 0x7F) >>  0);
						PRINT_LOG("MP3");
						*infoSize = (wk_long + 18);
					}
				}
			}
			
		} else if (infoFlacFile(filename, &inFmt, &inSample, &fileinfo) == TRUE) {
			PRINT_LOG("infoFlacFile");

			if (nChunk == 1) {
				memcpy(meta_list.h_list, "LIST", 4);
				meta_list.s_list = 0;
				memcpy(meta_list.h_info, "INFO", 4);
				memcpy(meta_list.h_isrc, "ISRC", 4);
				meta_list.s_isrc = 0;
				meta_list.p_isrc = NULL;
				memcpy(meta_list.h_icrd, "ICRD", 4);
				meta_list.s_icrd = 0;
				meta_list.p_icrd = NULL;
				memcpy(meta_list.h_iprd, "IPRD", 4);
				meta_list.s_iprd = 0;
				meta_list.p_iprd = NULL;
				memcpy(meta_list.h_ieng, "IENG", 4);
				meta_list.s_ieng = 0;
				meta_list.p_ieng = NULL;
				memcpy(meta_list.h_inam, "INAM", 4);
				meta_list.s_inam = 0;
				meta_list.p_inam = NULL;
				memcpy(meta_list.h_icop, "ICOP", 4);
				meta_list.s_icop = 0;
				meta_list.p_icop = NULL;
				memcpy(meta_list.h_ignr, "IGNR", 4);
				meta_list.s_ignr = 0;
				meta_list.p_ignr = NULL;
				memcpy(meta_list.h_isft, "ISFT", 4);
				meta_list.s_isft = 0;
				meta_list.p_isft = NULL;
				memcpy(meta_list.h_icmt, "ICMT", 4);
				meta_list.s_icmt = 0;
				meta_list.p_icmt = NULL;
				memcpy(meta_list.h_iart, "IART", 4);
				meta_list.s_iart = 0;
				meta_list.p_iart = NULL;

				fp = fopen(filename, "rb");
				if (fp == NULL) {
					break;
				}
				seekPtr = 4;
				for (end = 0; end == 0; ) {
					if (fseek_local(fp, seekPtr, SEEK_SET)) {
						break;
					}
					// ヘッダ読み込み
					rs = fread(header, 1, 4, fp);
					if (rs != 4) {
						break;
					}
					if ((header[0] & 0x80)) {
						end = 1;
					}
					wk_long = header[1] << 16 | header[2] << 8 | header[3];
					seekPtr += wk_long + 4;
					if ((header[0] & 0x7F) == 0x04) {
						PRINT_LOG("VORBIS_COMMENT");
						// VORBIS_COMMENT
						end = 1;
						// vendor_length
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = header[3] << 24 | header[2] << 16 | header[1]	<< 8 | header[0];
						if (fseek_local(fp, wk_long, SEEK_CUR)) {
							break;
						}
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						ovcom_count = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
						if (1) {
							char s[100];
							sprintf(s,"ovcom_count:%d",ovcom_count);
							PRINT_LOG(s);
						}
						for (i = 0; i < ovcom_count; i++) {
							rs = fread(header, 1, 4, fp);
							if (rs != 4) {
								break;
							}
							wk_long = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
							ovcom = (BYTE*)malloc(wk_long + 2);
							if (ovcom) {
								memset(ovcom, 0, wk_long + 2);
								rs = fread(ovcom, 1, wk_long, fp);
								if (wk_long != rs) {
									break;
								}
								p1 = (char*)ovcom;
								p2 = strchr((const char*)ovcom, '=');
								if (p2 != NULL) {
									*p2 = '\0';
									p2++;
								}
								if (stricmp(p1, "title") == 0 && strlen(p2) > 0) {
									if (meta_list.p_inam != NULL) {
										free(meta_list.p_inam);
										meta_list.p_inam = NULL;
										meta_list.s_inam = 0;
									}
									meta_list.p_inam = (BYTE*)malloc(strlen(p2));
									if (meta_list.p_inam != NULL) {
										meta_list.s_inam = strlen(p2);
										p_sjis = utf8_to_sjis(p2,&sjis_size);
										if (sjis_size > 0) {
											free(meta_list.p_inam);
											meta_list.p_inam = p_sjis;
											meta_list.s_inam = sjis_size;
										} else {
											memcpy(meta_list.p_inam, p2,meta_list.s_inam);
										}
									}
								} else if (stricmp(p1, "album") == 0 && strlen(p2) > 0) {
									if (meta_list.p_iprd != NULL) {
										free(meta_list.p_iprd);
										meta_list.p_iprd = NULL;
										meta_list.s_iprd = 0;
									}
									meta_list.p_iprd = (BYTE*)malloc(strlen(p2));
									if (meta_list.p_iprd != NULL) {
										meta_list.s_iprd = strlen(p2);
										p_sjis = utf8_to_sjis(p2,&sjis_size);
										if (sjis_size > 0) {
											free(meta_list.p_iprd);
											meta_list.p_iprd = p_sjis;
											meta_list.s_iprd = sjis_size;
										} else {
											memcpy(meta_list.p_iprd, p2,meta_list.s_iprd);
										}
									}
								} else if (stricmp(p1, "artist") == 0 && strlen(p2) > 0) {
									if (meta_list.p_iart != NULL) {
										free(meta_list.p_iart);
										meta_list.p_iart = NULL;
										meta_list.s_iart = 0;
									}
									meta_list.p_iart = (BYTE*)malloc(strlen(p2));
									if (meta_list.p_iart != NULL) {
										meta_list.s_iart = strlen(p2);
										if (1) {
											char s[100];
											sprintf(s,"IART:%s",p2);
											PRINT_LOG(s);
										}
										p_sjis = utf8_to_sjis(p2,&sjis_size);
										if (sjis_size > 0) {
if (1) {
	char s[100];
	sprintf(s,"IART(sjis):%s",p_sjis);
	PRINT_LOG(s);
}
											free(meta_list.p_iart);
											meta_list.p_iart = p_sjis;
											meta_list.s_iart = sjis_size;
										} else {
											memcpy(meta_list.p_iart, p2,meta_list.s_iart);
										}
									}
								} else if (stricmp(p1,"copyright") == 0 && strlen(p2) > 0) {
									if (meta_list.p_icop != NULL) {
										free(meta_list.p_icop);
										meta_list.p_icop = NULL;
										meta_list.s_icop = 0;
									}
									meta_list.p_icop = (BYTE*)malloc(strlen(p2));
									if (meta_list.p_icop != NULL) {
										meta_list.s_icop = strlen(p2);
										p_sjis = utf8_to_sjis(p2,&sjis_size);
										if (sjis_size > 0) {
											free(meta_list.p_icop);
											meta_list.p_icop = p_sjis;
											meta_list.s_icop = sjis_size;
										} else {
											memcpy(meta_list.p_icop, p2,meta_list.s_icop);
										}
									}

								} else if (stricmp(p1,"description") == 0 && strlen(p2) > 0) {
									if (meta_list.p_icmt != NULL) {
										free(meta_list.p_icmt);
										meta_list.p_icmt = NULL;
										meta_list.s_icmt = 0;
									}
									meta_list.p_icmt = (BYTE*)malloc(strlen(p2));
									if (meta_list.p_icmt != NULL) {
										meta_list.s_icmt = strlen(p2);
										p_sjis = utf8_to_sjis(p2,&sjis_size);
										if (p_sjis) {
											free(meta_list.p_icmt);
											meta_list.p_icmt = p_sjis;
											meta_list.s_icmt = sjis_size;
										} else {
											memcpy(meta_list.p_icmt, p2,meta_list.s_icmt);
										}
									}
								} else if (stricmp(p1, "date") == 0 && strlen(p2) > 0) {
									if (meta_list.p_icrd != NULL) {
										free(meta_list.p_icrd);
										meta_list.p_icrd = NULL;
										meta_list.s_icrd = 0;
									}
									meta_list.p_icrd = (BYTE*)malloc(strlen(p2));
									if (meta_list.p_icrd != NULL) {
										meta_list.s_icrd = strlen(p2);
										memcpy(meta_list.p_icrd, p2,meta_list.s_icrd);
									}
								}
								free(ovcom);
								ovcom = NULL;
							}
						}
						if (i == ovcom_count) {
							// サイズ計算
							wk_long = 0;
							if (meta_list.s_isrc > 0) {
								wk_long += meta_list.s_isrc + 8;
								if (meta_list.s_isrc & 0x01) wk_long++;
							}
							if (meta_list.s_icrd > 0) {
								wk_long += meta_list.s_icrd + 8;
								if (meta_list.s_icrd & 0x01) wk_long++;
							}
							if (meta_list.s_iprd > 0) {
								wk_long += meta_list.s_iprd + 8;
								if (meta_list.s_iprd & 0x01) wk_long++;
							}
							if (meta_list.s_ieng > 0) {
								wk_long += meta_list.s_ieng + 8;
								if (meta_list.s_ieng & 0x01) wk_long++;
							}
							if (meta_list.s_inam > 0) {
								wk_long += meta_list.s_inam + 8;
								if (meta_list.s_inam & 0x01) wk_long++;
							}
							if (meta_list.s_icop > 0) {
								wk_long += meta_list.s_icop + 8;
								if (meta_list.s_icop & 0x01) wk_long++;
							}
							if (meta_list.s_ignr > 0) {
								wk_long += meta_list.s_ignr + 8;
								if (meta_list.s_ignr & 0x01) wk_long++;
							}
							if (meta_list.s_isft > 0) {
								wk_long += meta_list.s_isft + 8;
								if (meta_list.s_isft & 0x01) wk_long++;
							}
							if (meta_list.s_icmt > 0) {
								wk_long += meta_list.s_icmt + 8;
								if (meta_list.s_icmt & 0x01) wk_long++;
							}
							if (meta_list.s_iart > 0) {
								wk_long += meta_list.s_iart + 8;
								if (meta_list.s_iart & 0x01) wk_long++;
							}
							meta_list.s_list = wk_long + 4;
							if (meta_list.s_list > 4) {
								*infoChunk = (char*)malloc(meta_list.s_list + 8 + 18);
								if (*infoChunk != NULL) {
									p1 = *infoChunk;
									memcpy(p1, meta_list.h_list, 12);
									p1 += 12;
									if (meta_list.s_isrc > 0) {
										memcpy(p1, meta_list.h_isrc, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_isrc,meta_list.s_isrc);
										p1 += meta_list.s_isrc;
										if (meta_list.s_isrc & 0x01) {
											*p1++ = '\0';
										}
									}
									if (meta_list.s_icrd > 0) {
										memcpy(p1, meta_list.h_icrd, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_icrd,meta_list.s_icrd);
										p1 += meta_list.s_icrd;
										if (meta_list.s_icrd & 0x01) {
											*p1++ = '\0';
										}
									}
									if (meta_list.s_iprd > 0) {
										memcpy(p1, meta_list.h_iprd, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_iprd,meta_list.s_iprd);
										p1 += meta_list.s_iprd;
										if (meta_list.s_iprd & 0x01) {
											*p1++ = '\0';
										}
									}
									if (meta_list.s_ieng > 0) {
										memcpy(p1, meta_list.h_ieng, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_ieng,meta_list.s_ieng);
										p1 += meta_list.s_ieng;
										if (meta_list.s_ieng & 0x01) {
											*p1++ = '\0';
										}
									}
									if (meta_list.s_inam > 0) {
										memcpy(p1, meta_list.h_inam, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_inam,meta_list.s_inam);
										p1 += meta_list.s_inam;
										if (meta_list.s_inam & 0x01) {
											*p1++ = '\0';
										}
									}
									if (meta_list.s_icop > 0) {
										memcpy(p1, meta_list.h_icop, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_icop,meta_list.s_icop);
										p1 += meta_list.s_icop;
										if (meta_list.s_icop & 0x01) {
											*p1++ = '\0';
										}
									}
									if (meta_list.s_ignr > 0) {
										memcpy(p1, meta_list.h_ignr, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_ignr,meta_list.s_ignr);
										p1 += meta_list.s_ignr;
										if (meta_list.s_ignr & 0x01) {
											*p1++ = '\0';
										}
									}
									if (meta_list.s_isft > 0) {
										memcpy(p1, meta_list.h_isft, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_isft,meta_list.s_isft);
										p1 += meta_list.s_isft;
										if (meta_list.s_isft & 0x01) {
											*p1++ = '\0';
										}
									}
									if (meta_list.s_icmt > 0) {
										memcpy(p1, meta_list.h_icmt, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_icmt,meta_list.s_icmt);
										p1 += meta_list.s_icmt;
										if (meta_list.s_icmt & 0x01) {
											*p1++ = '\0';
										}
									}
									if (meta_list.s_iart > 0) {
										memcpy(p1, meta_list.h_iart, 8);
										p1 += 8;
										memcpy(p1, meta_list.p_iart,meta_list.s_iart);
										p1 += meta_list.s_iart;
										if (meta_list.s_iart & 0x01) {
											*p1++ = '\0';
										}
									}
									*infoSize = ((p1 - *infoChunk));
									meta_list.s_list = *infoSize -8;
								}
							}
							break;
						}
						break;
					}
				}
			} else if (nChunk == 2) {
				memcpy(id3_d.talb_id,"TALB",4);		// アルバム名
				memcpy(id3_d.tpe1_id,"TPE1",4);		// 演奏者
				memcpy(id3_d.tpe2_id,"TPE2",4);		// バンド/オーケストラ
				memcpy(id3_d.comm_id,"COMM",4);		// コメント
				memcpy(id3_d.tpos_id,"TPOS",4);		// セット中の位置
				memcpy(id3_d.tcon_id,"TCON",4);		// 内容のタイプ
				memcpy(id3_d.tit2_id,"TIT2",4);		// 曲名
				memcpy(id3_d.trck_id,"TRCK",4);		// トラック番号
				memcpy(id3_d.tyer_id,"TYER",4);		// 年
				memcpy(id3_d.apic_id,"APIC",4);		// カバーアート

				// FlacファイルからID3タグの生成
				fp = fopen(filename, "rb");
				if (fp == NULL) {
					break;
				}
				seekPtr = 4;
				for (end = 0; end == 0; ) {
					if (fseek_local(fp, seekPtr, SEEK_SET)) {
						break;
					}
					// ヘッダ読み込み
					rs = fread(header, 1, 4, fp);
					if (rs != 4) {
						break;
					}
					if ((header[0] & 0x80)) {
						end = 1;
					}
					wk_long = (long)header[1] << 16 | (long)header[2] << 8 | (long)header[3];
					seekPtr += wk_long + 4;
					if ((header[0] & 0x7F) == 0x04) {
						// VORBIS_COMMENT

						// vendor_length
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = (long)header[3] << 24 | (long)header[2] << 16 | (long)header[1]	<< 8 | (long)header[0];
						if (fseek_local(fp, wk_long, SEEK_CUR)) {
							break;
						}
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						ovcom_count = (long)header[3] << 24 | (long)header[2] << 16 | (long)header[1] << 8 | (long)header[0];
						for (i = 0; i < ovcom_count; i++) {
							rs = fread(header, 1, 4, fp);
							if (rs != 4) {
								break;
							}
							wk_long = (long)header[3] << 24 | (long)header[2] << 16 | (long)header[1] << 8 | (long)header[0];
							ovcom = (BYTE*)malloc(wk_long + 2);
							if (ovcom) {
								memset(ovcom, 0, wk_long + 2);
								rs = fread(ovcom, 1, wk_long, fp);
								if (wk_long != rs) {
									break;
								}
								p1 = (char*)ovcom;
								p2 = strchr((const char*)ovcom, '=');
								if (p2 != NULL) {
									*p2 = '\0';
									p2++;
								}
								if (stricmp(p1, "title") == 0 && strlen(p2) > 0) {
									p_u16 = utf8_to_utf16(p2,&u16_size);
									if (p_u16) {
										id3_d.tit2_data = malloc(u16_size + 7);
										if (id3_d.tit2_data != NULL) {
											if (1) {
												char s[100];
												sprintf(s,"title:%d,%s",u16_size,p2);
												PRINT_LOG(s);
											}
											id3_d.tit2_data[0] = (BYTE)((u16_size + 1) >> 24);
											id3_d.tit2_data[1] = (BYTE)((u16_size + 1) >> 16);
											id3_d.tit2_data[2] = (BYTE)((u16_size + 1) >> 8);
											id3_d.tit2_data[3] = (BYTE)((u16_size + 1) >> 0);
											id3_d.tit2_data[4] = 0x00;
											id3_d.tit2_data[5] = 0x00;
											id3_d.tit2_data[6] = 0x01;
											memcpy(id3_d.tit2_data + 7,p_u16,u16_size);
											id3_d.tit2_size = u16_size + 7;
										}
										free(p_u16);
									}
								} else if (stricmp(p1, "album") == 0 && strlen(p2) > 0) {
									p_u16 = utf8_to_utf16(p2,&u16_size);
									if (p_u16) {
										id3_d.talb_data = malloc(u16_size + 7);
										if (id3_d.talb_data != NULL) {
											if (1) {
												char s[100];
												sprintf(s,"album:%d,%s",u16_size,p2);
												PRINT_LOG(s);
											}
											id3_d.talb_data[0] = (BYTE)((u16_size + 1) >> 24);
											id3_d.talb_data[1] = (BYTE)((u16_size + 1) >> 16);
											id3_d.talb_data[2] = (BYTE)((u16_size + 1) >> 8);
											id3_d.talb_data[3] = (BYTE)((u16_size + 1) >> 0);
											id3_d.talb_data[4] = 0x00;
											id3_d.talb_data[5] = 0x00;
											id3_d.talb_data[6] = 0x01;
											memcpy(id3_d.talb_data + 7,p_u16,u16_size);
											id3_d.talb_size = u16_size + 7;
										}
										free(p_u16);
									}
								} else if (stricmp(p1, "artist") == 0 && strlen(p2) > 0) {
									p_u16 = utf8_to_utf16(p2,&u16_size);
									if (p_u16) {
										id3_d.tpe1_data = malloc(u16_size + 7);
										if (id3_d.tpe1_data != NULL) {
											if (1) {
												char s[100];
												sprintf(s,"artist:%d,%s",u16_size,p2);
												PRINT_LOG(s);
											}
											id3_d.tpe1_data[0] = (BYTE)((u16_size + 1) >> 24);
											id3_d.tpe1_data[1] = (BYTE)((u16_size + 1) >> 16);
											id3_d.tpe1_data[2] = (BYTE)((u16_size + 1) >> 8);
											id3_d.tpe1_data[3] = (BYTE)((u16_size + 1) >> 0);
											id3_d.tpe1_data[4] = 0x00;
											id3_d.tpe1_data[5] = 0x00;
											id3_d.tpe1_data[6] = 0x01;
											memcpy(id3_d.tpe1_data + 7,p_u16,u16_size);
											id3_d.tpe1_size = u16_size + 7;
										}
										free(p_u16);
									}
								} else if (stricmp(p1,"description") == 0 && strlen(p2) > 0) {
									p_u16 = utf8_to_utf16(p2,&u16_size);
									if (p_u16) {
										id3_d.comm_data = malloc(u16_size + 7);
										if (id3_d.comm_data != NULL) {
											if (1) {
												char s[100];
												sprintf(s,"desc:%d,%s",u16_size,p2);
												PRINT_LOG(s);
											}
											id3_d.comm_data[0] = (BYTE)((u16_size + 1) >> 24);
											id3_d.comm_data[1] = (BYTE)((u16_size + 1) >> 16);
											id3_d.comm_data[2] = (BYTE)((u16_size + 1) >> 8);
											id3_d.comm_data[3] = (BYTE)((u16_size + 1) >> 0);
											id3_d.comm_data[4] = 0x00;
											id3_d.comm_data[5] = 0x00;
											id3_d.comm_data[6] = 0x01;
											memcpy(id3_d.comm_data + 7,p_u16,u16_size);
											id3_d.comm_size = u16_size + 7;
										}
										free(p_u16);
									}
								} else if (stricmp(p1, "date") == 0 && strlen(p2) > 0) {
									int s_len = strlen(p2) + 1;
									id3_d.tyer_data = malloc(5 + 7);
									if (id3_d.tyer_data != NULL) {
											if (1) {
												char s[100];
												sprintf(s,"date:%d,%s",s_len,p2);
												PRINT_LOG(s);
											}
										id3_d.tyer_data[0] = (BYTE)((s_len + 1) >> 24);
										id3_d.tyer_data[1] = (BYTE)((s_len + 1) >> 16);
										id3_d.tyer_data[2] = (BYTE)((s_len + 1) >> 8);
										id3_d.tyer_data[3] = (BYTE)((s_len + 1) >> 0);
										id3_d.tyer_data[4] = 0x00;
										id3_d.tyer_data[5] = 0x00;
										id3_d.tyer_data[6] = 0x00;
										id3_d.tyer_data[11] = 0x00;
										memcpy(id3_d.tyer_data + 7,p2,s_len);
										id3_d.tyer_size = s_len + 7;
									}
								} else if (stricmp(p1, "tracknumber") == 0 && strlen(p2) > 0) {
									int s_len = strlen(p2) + 1;
									id3_d.trck_data = malloc(s_len + 7);
									if (id3_d.trck_data != NULL) {
											if (1) {
												char s[100];
												sprintf(s,"tracknumber:%d,%s",s_len,p2);
												PRINT_LOG(s);
											}
										id3_d.trck_data[0] = (BYTE)((s_len + 1) >> 24);
										id3_d.trck_data[1] = (BYTE)((s_len + 1) >> 16);
										id3_d.trck_data[2] = (BYTE)((s_len + 1) >> 8);
										id3_d.trck_data[3] = (BYTE)((s_len + 1) >> 0);
										id3_d.trck_data[4] = 0x00;
										id3_d.trck_data[5] = 0x00;
										id3_d.trck_data[6] = 0x00;
										memcpy(id3_d.trck_data + 7,p2,s_len);
										id3_d.trck_size = s_len + 7;
									}
								}
								free(ovcom);
								ovcom = NULL;
							}
						}
					} else if ((header[0] & 0x7F) == 0x06) {
						// PICTURE
						PRINT_LOG("PICTURE");
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = (long)header[0] << 24 | (long)header[1] << 16 | (long)header[2]	<< 8 | (long)header[3];
						flac_pic.type = wk_long;
						PRINT_LOG("pic.type");

						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = (long)header[0] << 24 | (long)header[1] << 16 | (long)header[2]	<< 8 | (long)header[3];
						flac_pic.mime_size = wk_long;
						PRINT_LOG("mime_size");
						if (flac_pic.mime_size > 0) {
							flac_pic.mime_data = malloc(flac_pic.mime_size + 1);
							if (flac_pic.mime_data != NULL) {
								memset(flac_pic.mime_data,0,flac_pic.mime_size + 1);
								rs = fread(flac_pic.mime_data, 1, flac_pic.mime_size, fp);
								if (rs != flac_pic.mime_size) {
									break;
								}
							} else {
								break;
							}
						}
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = (long)header[0] << 24 | (long)header[1] << 16 | (long)header[2]	<< 8 | (long)header[3];
						flac_pic.desc_size = wk_long;
						PRINT_LOG("desc_size");
						if (flac_pic.desc_size > 0) {
							flac_pic.desc_data = malloc(flac_pic.desc_size + 1);
							if (flac_pic.desc_data != NULL) {
								memset(flac_pic.desc_data,0,flac_pic.desc_size + 1);
								rs = fread(flac_pic.desc_data, 1, flac_pic.desc_size, fp);
								if (rs != flac_pic.desc_size) {
									break;
								}
								flac_pic.desc_size = flac_pic.desc_size + 1;
							} else {
								break;
							}
						} else {
							flac_pic.desc_data = malloc(1);
							if (flac_pic.desc_data != NULL) {
								flac_pic.desc_data[0] = 0x00;
								flac_pic.desc_size = 1;
							}
						}
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = (long)header[0] << 24 | (long)header[1] << 16 | (long)header[2]	<< 8 | (long)header[3];
						flac_pic.width = wk_long;
						PRINT_LOG("width");

						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = (long)header[0] << 24 | (long)header[1] << 16 | (long)header[2]	<< 8 | (long)header[3];
						flac_pic.height = wk_long;
						PRINT_LOG("height");

						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = (long)header[0] << 24 | (long)header[1] << 16 | (long)header[2]	<< 8 | (long)header[3];
						flac_pic.depth  = wk_long;
						PRINT_LOG("depth");

						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = (long)header[0] << 24 | (long)header[1] << 16 | (long)header[2]	<< 8 | (long)header[3];
						flac_pic.color  = wk_long;
						PRINT_LOG("color");

						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = (long)header[0] << 24 | (long)header[1] << 16 | (long)header[2]	<< 8 | (long)header[3];
						flac_pic.pic_data_size = wk_long;
						PRINT_LOG("data_size");
						if (flac_pic.pic_data_size > 0) {
							flac_pic.pic_data = malloc(flac_pic.pic_data_size);
							if (flac_pic.pic_data != NULL) {
								rs = fread(flac_pic.pic_data, 1, flac_pic.pic_data_size,fp);
								if (rs != flac_pic.pic_data_size) {
									break;
								}
							} else {
								break;
							}
						}
						end = 1;
						if (1) {
							char s[2048];
							sprintf(s,"FlacPic:type:%d,mime:%d,%s,desc:%d,width:%d,height:%d,depth:%d,color:%d,pic_data_size:%d",flac_pic.type,flac_pic.mime_size,flac_pic.mime_data,flac_pic.desc_size,flac_pic.width,flac_pic.height,flac_pic.depth,flac_pic.color,flac_pic.pic_data_size);
							PRINT_LOG(s);
						}
						if (flac_pic.mime_size > 0 && flac_pic.width > 0 && flac_pic.height > 0 && flac_pic.pic_data_size > 0) {
							PRINT_LOG("apic set");
							wk_long = 6;	// size + flag(2)
							wk_long += (1 + flac_pic.mime_size + 1);
							wk_long += 1;	// picture type
							wk_long += flac_pic.desc_size;
							wk_long += flac_pic.pic_data_size;
							id3_d.apic_data = malloc(wk_long + 1);
							if (id3_d.apic_data != NULL) {
								id3_d.apic_data[0] = (BYTE)((wk_long - 6) >> 24);
								id3_d.apic_data[1] = (BYTE)((wk_long - 6) >> 16);
								id3_d.apic_data[2] = (BYTE)((wk_long - 6) >> 8);
								id3_d.apic_data[3] = (BYTE)((wk_long - 6) >> 0);
								id3_d.apic_data[4] = 0x00;	// Flag
								id3_d.apic_data[5] = 0x00;	// Flag
								id3_d.apic_data[6] = 0x00;	// Text Encording
								p1 = &id3_d.apic_data[7];
								strcpy(p1,flac_pic.mime_data);p1 += (strlen(flac_pic.mime_data) + 1);
								*p1 = (BYTE)flac_pic.type;p1++;
								memcpy(p1,flac_pic.desc_data,flac_pic.desc_size);p1 += flac_pic.desc_size;
								memcpy(p1,flac_pic.pic_data,flac_pic.pic_data_size);
								id3_d.apic_size = wk_long;
							} else {
								break;
							}
						}
					}
				}
				wk_long = 0;
				if (id3_d.talb_size > 0) {wk_long += 4;wk_long += id3_d.talb_size;}
				if (id3_d.tpe1_size > 0) {wk_long += 4;wk_long += id3_d.tpe1_size;}
				if (id3_d.tpe2_size > 0) {wk_long += 4;wk_long += id3_d.tpe2_size;}
				if (id3_d.comm_size > 0) {wk_long += 4;wk_long += id3_d.comm_size;}
				if (id3_d.tpos_size > 0) {wk_long += 4;wk_long += id3_d.tpos_size;}
				if (id3_d.tcon_size > 0) {wk_long += 4;wk_long += id3_d.tcon_size;}
				if (id3_d.tit2_size > 0) {wk_long += 4;wk_long += id3_d.tit2_size;}
				if (id3_d.trck_size > 0) {wk_long += 4;wk_long += id3_d.trck_size;}
				if (id3_d.tyer_size > 0) {wk_long += 4;wk_long += id3_d.tyer_size;}
				if (id3_d.apic_size > 0) {wk_long += 4;wk_long += id3_d.apic_size;}
				if (1) {
					char s[100];
					sprintf(s,"ID3:%d",wk_long);
					PRINT_LOG(s);
				}
				if (wk_long > 0) {
					id3_h.id3_tag[0] = 'I';
					id3_h.id3_tag[1] = 'D';
					id3_h.id3_tag[2] = '3';
					id3_h.id3_ver[0] = 3;
					b3 = (BYTE)(wk_long);b3 &= 0x7F;
					b2 = (BYTE)(wk_long >> 7);b2 &= 0x7F;
					b1 = (BYTE)(wk_long >> 14);b1 &= 0x7F;
					b0 = (BYTE)(wk_long >> 21);b0 &= 0x7F;
					id3_h.id3_size[0] = b0;
					id3_h.id3_size[1] = b1;
					id3_h.id3_size[2] = b2;
					id3_h.id3_size[3] = b3;

					*infoChunk = (char*)malloc(wk_long + 18);
					if (*infoChunk != NULL) {
						p1 = *infoChunk + 8;
						memcpy(p1,&id3_h,10);p1 += 10;
						if (id3_d.talb_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"talb:%d",id3_d.talb_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.talb_id,4);p1 += 4;
							memcpy(p1,id3_d.talb_data,id3_d.talb_size);
							p1 += id3_d.talb_size;
						}
						if (id3_d.tpe1_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"tpe1:%d",id3_d.tpe1_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.tpe1_id,4);p1 += 4;
							memcpy(p1,id3_d.tpe1_data,id3_d.tpe1_size);
							p1 += id3_d.tpe1_size;
						}
						if (id3_d.tpe2_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"tpe2:%d",id3_d.tpe2_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.tpe2_id,4);p1 += 4;
							memcpy(p1,id3_d.tpe2_data,id3_d.tpe2_size);
							p1 += id3_d.tpe2_size;
						}
						if (id3_d.comm_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"comm:%d",id3_d.comm_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.comm_id,4);p1 += 4;
							memcpy(p1,id3_d.comm_data,id3_d.comm_size);
							p1 += id3_d.comm_size;
						}
						if (id3_d.tpos_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"tpos:%d",id3_d.tpos_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.tpos_id,4);p1 += 4;
							memcpy(p1,id3_d.tpos_data,id3_d.tpos_size);
							p1 += id3_d.tpos_size;
						}
						if (id3_d.tcon_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"tcon:%d",id3_d.tcon_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.tcon_id,4);p1 += 4;
							memcpy(p1,id3_d.tcon_data,id3_d.tcon_size);
							p1 += id3_d.tcon_size;
						}
						if (id3_d.tit2_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"tit2:%d",id3_d.tit2_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.tit2_id,4);p1 += 4;
							memcpy(p1,id3_d.tit2_data,id3_d.tit2_size);
							p1 += id3_d.tit2_size;
						}
						if (id3_d.trck_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"trck:%d",id3_d.trck_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.trck_id,4);p1 += 4;
							memcpy(p1,id3_d.trck_data,id3_d.trck_size);
							p1 += id3_d.trck_size;
						}
						if (id3_d.tyer_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"tyer:%d",id3_d.tyer_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.tyer_id,4);p1 += 4;
							memcpy(p1,id3_d.tyer_data,id3_d.tyer_size);
							p1 += id3_d.tyer_size;
						}
						if (id3_d.apic_size > 0) {
							if (1) {
								char s[100];
								sprintf(s,"apic:%d",id3_d.apic_size);
								PRINT_LOG(s);
							}
							memcpy(p1,id3_d.apic_id,4);p1 += 4;
							memcpy(p1,id3_d.apic_data,id3_d.apic_size);
							p1 += id3_d.apic_size;
						}
						*infoSize = (p1 - *infoChunk);
						p1 = *infoChunk;
						p1[0] = 'i';
						p1[1] = 'd';
						p1[2] = '3';
						p1[3] = ' ';
						p1[4] = (BYTE)((*infoSize - 8));
						p1[5] = (BYTE)((*infoSize - 8) >>  8);
						p1[6] = (BYTE)((*infoSize - 8) >> 16);
						p1[7] = (BYTE)((*infoSize - 8) >> 24);
						if (*infoSize & 0x01) *infoSize = *infoSize + 1;
						if (1) {
							char s[200];
							sprintf(s,"ID3 size:%ld",*infoSize);
							PRINT_LOG(s);
						}
					}
				}
			} else {
				break;
			}
		}
	}
	while (0);

	if (ovcom != NULL) {
		free(ovcom);
	}

	if (meta_list.p_isrc != NULL)
		free(meta_list.p_isrc);
	if (meta_list.p_icrd != NULL)
		free(meta_list.p_icrd);
	if (meta_list.p_iprd != NULL)
		free(meta_list.p_iprd);
	if (meta_list.p_ieng != NULL)
		free(meta_list.p_ieng);
	if (meta_list.p_inam != NULL)
		free(meta_list.p_inam);
	if (meta_list.p_icop != NULL)
		free(meta_list.p_icop);
	if (meta_list.p_ignr != NULL)
		free(meta_list.p_ignr);
	if (meta_list.p_isft != NULL)
		free(meta_list.p_isft);
	if (meta_list.p_icmt != NULL)
		free(meta_list.p_icmt);
	if (meta_list.p_iart != NULL)
		free(meta_list.p_iart);

	if (flac_pic.mime_data != NULL) free(flac_pic.mime_data);
	if (flac_pic.pic_data  != NULL) free(flac_pic.pic_data);

	if (id3_d.talb_data != NULL) free(id3_d.talb_data);
	if (id3_d.tpe1_data != NULL) free(id3_d.tpe1_data);
	if (id3_d.tpe2_data != NULL) free(id3_d.tpe2_data);
	if (id3_d.comm_data != NULL) free(id3_d.comm_data);
	if (id3_d.tpos_data != NULL) free(id3_d.tpos_data);
	if (id3_d.tcon_data != NULL) free(id3_d.tcon_data);
	if (id3_d.tit2_data != NULL) free(id3_d.tit2_data);
	if (id3_d.trck_data != NULL) free(id3_d.trck_data);
	if (id3_d.apic_data != NULL) free(id3_d.apic_data);

	if (fp != NULL) {
		fclose(fp);
	}

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_GetChunkInfo
 * Description : Wave ファイルのChunk情報を返す
 * ---
 *	filename   : ファイル名
 *	info	   : 情報を格納するバッファ
 *	infoSize   : info のバイト数
 *
 */
int PLG_GetChunkInfo(char *filename, char *info, int infoSize)
/*
 * チャンク情報の取得
 */
{
	FILE *fp;
	WFEX wfex;
	BROADCAST_EXT *bext;
	SOUNDFMT inFmt;
	DWORD inSample;
	int retValue;
	long chunkSize;
	DWORD rs;
	DWORD dw_temp;
	int remInfoSize;
	char work[300];
	char *listData;
	long listSize;
	char *data;
	char *p, *cr;
	long len;
	int i;
	__int64 offset;
	__int64 timeRef;
	int key_writed;
	__int64 seekPtr;
	int end;
	BYTE header[40];
	long wk_long;
	long ovcom_count;
	BYTE *ovcom;
	char *p1, *p2;
	DSF dsf;
	DSF_FMT dsf_fmt;
	DSF_DATA dsf_data;
	RF64_DS64 ds64;
	int index;
	char tag[5];
	int flg_wave_file;
	int flg_rf64_file;

	retValue = STATUS_FILE_READ_ERR;
	data = NULL;
	listData = NULL;
	bext = NULL;
	ovcom = NULL;
	end = 0;
	do {
		remInfoSize = infoSize;
		memset(info, 0, infoSize);

		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		retValue = STATUS_UNKNOWN_FORMAT;

		flg_wave_file = flg_rf64_file = FALSE;
		flg_wave_file = infoWavFile(filename, &inFmt, &inSample, NULL);
		if (flg_wave_file == FALSE) {
			flg_rf64_file = infoRF64File(filename, &inFmt, &inSample, NULL);
		}
		if (flg_wave_file) {
			// Wave File
			strcpy(info, "AWAVE\n"); // Add "Wave"
			sprintf(work, "IOffset\t-\n");
			strcat(info, work);
			sprintf(work, "ISize\t-\n");
			sprintf(work, "D\n");	// Down Tree
			strcat(info, work);
		}
		if (flg_rf64_file) {
			// Wave File
			strcpy(info, "ARF64\n"); // Add "RF64"
			sprintf(work, "IOffset\t-\n");
			strcat(info, work);
			sprintf(work, "ISize\t-\n");
			sprintf(work, "D\n");	// Down Tree
			strcat(info, work);
		}
		if (flg_wave_file || flg_rf64_file) {
			index = 1;
			tag[4] = '\0';
			while (findWavChunk_N(fp,index,tag,&chunkSize)) {
				if (strcmp("fmt ",tag) == 0) {
					strcat(info, "Afmt\n"); // Add "fmt"
					offset = ftell_local(fp);
					rs = fread(&wfex, 1, chunkSize, fp);
					if (rs != chunkSize) {
						break;
					}
					sprintf(work, "IID\tfmt \n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", rs + 8);
					strcat(info, work);

					sprintf(work, "IwFormatTag\t%04x", wfex.Format.wFormatTag);
					strcat(info, work);
					if (wfex.Format.wFormatTag == WF_PCM) {
						strcat(info, "(WAVE_FORMAT_PCM)\n");
					} else if (wfex.Format.wFormatTag == WF_IEEE_FLOAT) {
						strcat(info, "(WAVE_FORMAT_IEEE_FLOAT)\n");
					} else if (wfex.Format.wFormatTag == 0xFFFE) {
						strcat(info, "(WAVE_FORMAT_EXTENSIBLE)\n");
					} else {
						strcat(info, "\n");
					}
					sprintf(work, "InChannels\t%d\n", wfex.Format.nChannels);
					strcat(info, work);
					sprintf(work, "InSamplesPerSec\t%ld\n",wfex.Format.nSamplesPerSec);
					strcat(info, work);
					sprintf(work, "IwBitsPerSample\t%d\n",wfex.Format.wBitsPerSample);
					strcat(info, work);
					sprintf(work, "InBlockAlign\t%d\n", wfex.Format.nBlockAlign);
					strcat(info, work);
					sprintf(work, "InAvgBytesPerSec\t%ld\n",wfex.Format.nAvgBytesPerSec);
					strcat(info, work);
					if (wfex.Format.wFormatTag == 0xFFFE) {
						sprintf(work, "IwValidBitsPerSample\t%d\n",wfex.Samples.wValidBitsPerSample);
						strcat(info, work);
						sprintf(work, "IwSamplesPerBlock\t%d\n",wfex.Samples.wSamplesPerBlock);
						strcat(info, work);
						sprintf(work, "IdwChannelMask\t%X\n", wfex.dwChannelMask);
						strcat(info, work);
					}
				}
				if (strcmp(tag,"data") == 0) {
					strcat(info, "Adata\n"); // Add "data"
					offset = ftell_local(fp);
					sprintf(work, "IID\tdata\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", chunkSize + 8);
					strcat(info, work);

					if (chunkSize > 0L) {
						dw_temp = chunkSize;
						dw_temp /= (wfex.Format.wBitsPerSample / 8);
						dw_temp /= (wfex.Format.nChannels);
						sprintf(work, "ISampleCount\t%ld\n", dw_temp);
						strcat(info, work);
					}
				}
				if (strcmp(tag,"LIST") == 0) {
					if (chunkSize > 20) {
						offset = ftell_local(fp); listSize = chunkSize;
						listData = malloc(chunkSize);
						if (listData != NULL) {
							rs = fread(listData, 1, listSize, fp);
							if (rs != listSize) {
								free(listData);
								listData = NULL;
							} else {
								if (memcmp(listData, "INFO", 4) != 0) {
									free(listData);
									listData = NULL;
								}
							}
						}
					}
					if (listData != NULL) {
						strcat(info, "ALIST(INFO)\n");
						// Add "INFO"
						sprintf(work, "IID\tLIST(INFO)\n");
						strcat(info, work);
						sprintf(work, "IOffset\t%08lXh\n", offset - 8);
						strcat(info, work);
						sprintf(work, "ISize\t%llu\n", listSize + 8);
						strcat(info, work); p = listData + 4; listSize -= 4;
						while (listSize >= 8) {
							strcat(info, "I");
							memset(work, 0, 5);
							strncpy(work, p, 4);
							strcat(info, work);
							strcat(info, "\t");
							p += 4;
							len = *(long*)p;
							p += 4;
							memset(work, 0, 256);
							strncpy(work, p, len < 256 ? len : 255);
							if (len & 1) {
								len++;
							}
							strcat(info, work);
							strcat(info, "\n");
							p += len;
							if (listSize <= len + 8) {
								break;
							}
							listSize -= (len + 8);
						}
						free(listData);
						listData = NULL;
					}
				}
				if (strcmp(tag,"JUNK") == 0) {
					strcat(info, "AJUNK\n"); // Add "data"
					offset = ftell_local(fp);
					sprintf(work, "IID\tJUNK\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%ld\n", chunkSize + 8);
					strcat(info, work);
				}
				if (strcmp(tag,"id3 ") == 0) {
					strcat(info, "Aid3 \n"); // Add "data"
					offset = ftell_local(fp);
					sprintf(work, "IID\tID3\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", chunkSize + 8);
					strcat(info, work);
				}
				if (strcmp(tag,"bext") == 0) {
					// BWF
					strcat(info, "Abext\n"); // Add "bext"
					offset = ftell_local(fp);
					bext = (BROADCAST_EXT*)malloc(chunkSize);
					if (bext != NULL) {
						rs = fread(bext, 1, chunkSize, fp);
						if (rs != chunkSize) {
							free(bext); 
							bext = NULL;
						}
					}
					if (bext != NULL) {
						sprintf(work, "IID\tbext\n");
						strcat(info, work);
						sprintf(work, "IOffset\t%08lXh\n", offset - 8);
						strcat(info, work);
						sprintf(work, "ISize\t%llu\n", chunkSize + 8);
						strcat(info, work); strcat(info, "IDescription\t");
						strcpy(work, "-");
						if (strlen(bext->description) > 0) {
							memset(work, 0, 257);
							strncpy(work, bext->description, 256);
						}
						strcat(info, work);
						strcat(info, "\n");
						strcat(info, "IOriginator\t");
						strcpy(work, "-");
						if (strlen(bext->originator) > 0) {
							memset(work, 0, 257);
							strncpy(work, bext->originator, 32);
						}
						strcat(info, work); strcat(info, "\n");

						strcat(info, "IOriginatorReference\t");
						strcpy(work, "-");
						if (strlen(bext->originatorReference) > 0) {
							memset(work, 0, 257);
							strncpy(work, bext->originatorReference, 32);
						}
						strcat(info, work); 
						strcat(info, "\n");

						strcat(info, "IOriginationDate\t"); 
						strcpy(work, "-");
						if (strlen(bext->originationDate) > 0) {
							memset(work, 0, 257);
							strncpy(work, bext->originationDate, 10);
						}
						strcat(info, work);
						strcat(info, "\n");

						strcat(info, "IOriginationTime\t");
						strcpy(work, "-");
						if (strlen(bext->originationTime) > 0) {
							memset(work, 0, 257);
							strncpy(work, bext->originationTime, 8);
						}
						strcat(info, work); 
						strcat(info, "\n");

						strcat(info, "ITimeReference\t");
						timeRef = bext->timeReferenceHigh; timeRef <<= 32;
						timeRef |= bext->timeReferenceLow;
						timeRef /= wfex.Format.nSamplesPerSec;
						sprintf(work, "%d:%02d:%02d", (int)(timeRef / 3600),(int)((timeRef / 60) % 60), (int)(timeRef % 60));
						strcat(info, work); strcat(info, "\n");

						strcat(info, "IVersion\t");
						sprintf(work, "%d", bext->version);
						strcat(info, work); strcat(info, "\n");

						strcat(info, "IUMID\t"); strcpy(work, "-");
						if (strlen((char*)bext->UMID) > 0) {
							memset(work, 0, 257);
							strncpy(work, (char*)bext->UMID, 64);
						}
						strcat(info, work); strcat(info, "\n");

						strcat(info, "ICodingHistory\t"); strcpy(work, "-");
						key_writed = 1;
						if (strlen(bext->codingHistory) > 0) {
							p = bext->codingHistory;
							while (strlen(p) > 0) {
								if (key_writed == 0) {
									strcat(info, "\nI \t");
								}
								cr = strchr(p, '\n');
								if (cr) {
									*cr = '\0';
								}
								memset(work, 0, 257);
								strncpy(work, p, 256);
								strcat(info, work);
								if (cr) {
									p = (cr + 1);
								}
								key_writed = 0;
							}
						}
						if (key_writed == 1) {
							strcat(info, work);
						}
						strcat(info, "\n"); 
						free(bext); 
						bext = NULL;
					}
				}
				if (strcmp(tag,"qlty") == 0) {
					// BWF(qlty)
					strcat(info, "Aqlty\n"); // Add "qlty"
					offset = ftell_local(fp);
					sprintf(work, "IID\tqlty\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", chunkSize + 8);
					strcat(info, work);
				}
				if (strcmp(tag,"levl") == 0) {
					// BWF(qlty)
					strcat(info, "Alevl\n"); // Add "levl"
					offset = ftell_local(fp);
					sprintf(work, "IID\tlevl\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", chunkSize + 8);
					strcat(info, work);
				}
				if (strcmp(tag,"link") == 0) {
					// BWF(link)
					strcat(info, "Alink\n"); // Add "link"
					offset = ftell_local(fp);
					sprintf(work, "IID\tlink\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", chunkSize + 8);
					strcat(info, work);
				}
				if (strcmp(tag,"axml") == 0) {
					// BWF(axml)
					strcat(info, "Aaxml\n"); // Add "axml"
					offset = ftell_local(fp);
					sprintf(work, "IID\taxml\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", chunkSize + 8);
					strcat(info, work);
				}
				if (strcmp(tag,"cue ") == 0) {
					// BWFJ(cue)
					strcat(info, "Acue\n"); // Add "cue"
					offset = ftell_local(fp);
					sprintf(work, "IID\tcue \n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", chunkSize + 8);
					strcat(info, work);
				}
				if (strcmp(tag,"plst") == 0) {
					// BWFJ(plst)
					strcat(info, "Aplst\n"); // Add "plst"
					offset = ftell_local(fp);
					sprintf(work, "IID\tplst\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", chunkSize + 8);
					strcat(info, work);
				}
				if (strcmp(tag,"adtl") == 0) {
					// BWFJ Associated data chunk(adtl)
					strcat(info, "ALIST(adtl)\n"); // Add "adtl"
					offset = ftell_local(fp);
					sprintf(work, "IID\tadtl\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%llu\n", chunkSize + 8);
					strcat(info, work);
				}
				if (strcmp(tag,"ds64") == 0) {
					__int64 qw; 
					strcat(info, "Ads64\n"); // Add "fmt"
					offset = ftell_local(fp);
					if (fseek_local(fp,-8,SEEK_CUR)) {
						break;
					}
					rs = fread(&ds64, 1, chunkSize + 8, fp);
					if (rs != chunkSize + 8) {
						break;
					}
					sprintf(work, "IID\tds64\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset);
					strcat(info, work); sprintf(work, "ISize\t%llu\n", rs);
					strcat(info, work); qw = ds64.riffSizeHigh; qw <<= 32;
					qw |= ds64.riffSizeLow;
#ifdef __GNUC__
					sprintf(work, "IriffSize\t%llu\n", qw);
#else
					sprintf(work, "IriffSize\t%I64d\n", qw);
#endif
					strcat(info, work);

					qw = ds64.dataSizeHigh; qw <<= 32; qw |= ds64.dataSizeLow;
#ifdef __GNUC__
					sprintf(work, "IdataSize\t%llu\n", qw);
#else
					sprintf(work, "IdataSize\t%I64d\n", qw);
#endif
					strcat(info, work);

					qw = ds64.sampleCountHigh; qw <<= 32;
					qw |= ds64.sampleCountLow;
#ifdef __GNUC__
					sprintf(work, "IsampleCount\t%llu\n", qw);
#else
					sprintf(work, "IsampleCount\t%I64d\n", qw);
#endif
					strcat(info, work);
					sprintf(work, "ItableLength\t%ld\n", ds64.tableLength);
					strcat(info, work);
				}
				index++;
			}
			retValue = STATUS_SUCCESS;

			// List(info)

		}
		if (infoFlacFile(filename, &inFmt, &inSample, NULL)) {
			strcpy(info, "AFLAC\n"); // Add "Flac"
			sprintf(work, "IOffset\t-\n");
			strcat(info, work);
			sprintf(work, "ISize\t-\n");
			strcat(info, work);
			seekPtr = 4;
			if (fseek_local(fp, seekPtr, SEEK_SET)) {
				break;
			}
			rs = fread(header, 1, 38, fp);
			if (rs != 38) {
				break;
			}
			if ((header[0] & 0x80) != 0x00) {
				end = 1;
			}
			if ((header[0] & 0x7F) == 0x00) {
				//
				long sr;
				long ch;
				long bit;
				DWORD dw;
				long ld;
				__int64 total;
				strcat(info, "D\n"); // Down Tree
				strcat(info, "ASTREAMINFO\n"); // Add "STREAMINFO"
				offset = 4;
				sprintf(work, "IBlock type\t%d\n", header[0] & 0x7F);
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset);
				strcat(info, work);
				ld = header[1];
				ld <<= 8;ld |= header[2];
				ld <<= 8;ld |= header[3];
				sprintf(work, "ISize\t%ld\n", ld + 4); strcat(info, work);

				ld = header[4];
				ld <<= 8;ld |= header[5];

				sprintf(work, "Iminimum block size\t%ld\n", ld);
				strcat(info, work);

				ld = header[6]; ld <<= 8; ld |= header[7];

				sprintf(work, "Imaximum block size\t%ld\n", ld);
				strcat(info, work);

				ld = header[8]; ld <<= 8; ld |= header[9]; ld <<= 8;
				ld |= header[10];

				sprintf(work, "Iminimum frame size\t%ld\n", ld);
				strcat(info, work);

				ld = header[11]; ld <<= 8; ld |= header[12]; ld <<= 8;
				ld |= header[13];

				sprintf(work, "Imaximum frame size\t%ld\n", ld);
				strcat(info, work);

				sr = header[14]; sr <<= 8; sr |= header[15]; sr <<= 4;
				sr |= (header[16] >> 4) & 0x0F;

				sprintf(work, "Isample rate\t%ld\n", sr);
				strcat(info, work);

				ch = (header[16] >> 1) & 0x07; ch++;

				sprintf(work, "Ichannels\t%ld\n", ch); strcat(info, work);

				bit = (header[16] & 0x01) << 4;
				bit |= (header[17] >> 4) & 0x0F; bit++;

				sprintf(work, "Ibits per sample\t%ld\n", bit);
				strcat(info, work);

				total = (header[17] & 0x0F); total <<= 8;
				total |= header[18]; total <<= 8; total |= header[19];
				total <<= 8; total |= header[20];

#ifdef __GUNC__
				sprintf(work, "Itotal samples\t%lld\n", total);
#else
				sprintf(work, "Itotal samples\t%I64d\n", total);
#endif
				strcat(info, work);

				strcat(info, "IMD5\t"); 
				memset(work, 0, 300); 
				p1 = work;
				for (i = 0; i < 16; i++) {
					sprintf(p1, "%02x", header[21 + i]);
					p1 = work + strlen(work);
				}
				strcat(p1, "\n"); 
				strcat(info, work);
				retValue = STATUS_SUCCESS;
			}

			seekPtr = 4;
			for (; end == 0; ) {
				if (fseek_local(fp, seekPtr, SEEK_SET)) {
					break;
				}
				// ヘッダ読み込み
				rs = fread(header, 1, 4, fp);
				if (rs != 4) {
					break;
				}
				if ((header[0] & 0x80)) {
					end = 1;
				}
				wk_long = header[1] << 16 | header[2] << 8 | header[3];
				offset = seekPtr; seekPtr += wk_long + 4;
				if (!((header[0] & 0x7F) == 0x00 || (header[0] & 0x7F) == 0x04)) {
					long ld;
					wk_long = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
					if ((header[0] & 0x7F) == 1) {
						strcat(info, "APADDING\n"); // Add "PADDING"
					} else if ((header[0] & 0x7F) == 2) {
						strcat(info, "AAPPLICATION\n"); // Add "APPLICATION"
					} else if ((header[0] & 0x7F) == 3) {
						strcat(info, "ASEEKTABLE\n"); // Add "APPLICATION"
					} else if ((header[0] & 0x7F) == 5) {
						strcat(info, "ACUESHEET\n"); // Add "CUESHEET"
					} else if ((header[0] & 0x7F) == 6) {
						strcat(info, "APICTURE\n"); // Add "PICTURE"
					}
					strcat(info, work);
					sprintf(work, "IBlock type\t%d\n", header[0] & 0x7F);
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset);
					strcat(info, work);
					sprintf(work, "ISize\t%ld\n", wk_long);
					strcat(info, work);
				}
				if ((header[0] & 0x7F) == 0x04) {
					// VORBIS_COMMENT
					// end = 1;
					// vendor_length
					rs = fread(header, 1, 4, fp);
					if (rs != 4) {
						break;
					}
					wk_long = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
					if (fseek_local(fp, wk_long, SEEK_CUR)) {
						break;
					}
					rs = fread(header, 1, 4, fp);
					if (rs != 4) {
						break;
					}
					ovcom_count = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
					if (ovcom_count > 0){
						long ld;
						strcat(info, "AVORBIS_COMMENT\n");
						// Add "VORBIS_COMMENT"
						strcat(info, work);
						sprintf(work, "IBlock type\t%d\n",header[0] & 0x7F); strcat(info, work);
						sprintf(work, "IOffset\t%08lXh\n", offset);
						strcat(info, work);
						sprintf(work, "ISize\t%ld\n", wk_long);
						strcat(info, work);
					}

					for (i = 0; i < ovcom_count; i++) {
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
						ovcom = (BYTE*)malloc(wk_long + 2);
						if (ovcom) {
							memset(ovcom, 0, wk_long + 2);
							rs = fread(ovcom, 1, wk_long, fp);
							if (wk_long != rs) {
								break;
							}
							p1 = (char*)ovcom;
							p2 = strchr((const char*)ovcom, '=');
							if (p2 != NULL) {
								*p2 = '\0';
								p2++;
							}
							if (p1 != NULL && p2 != NULL && strlen(p1) > 0 && strlen(p2) > 0) {
								sprintf(work, "I%s\t%s\n", p1, p2);
								strcat(info, work);
							}
							free(ovcom); 
							ovcom = NULL;
						}
					}
				}
			}
		}
		if (infoDsfFile(filename, &inFmt, &inSample, NULL)) {
			strcpy(info, "ADSF\n"); // Add "DSF"
			sprintf(work, "IOffset\t-\n"); strcat(info, work);
			sprintf(work, "ISize\t-\n"); strcat(info, work);

			seekPtr = 0;
			if (fseek_local(fp, seekPtr, SEEK_SET)) {
				break;
			}
			// DSF
			rs = fread(&dsf, 1, sizeof(DSF), fp);
			if (rs != sizeof(DSF)) {
				break;
			}
			strcat(info, "D\n"); // Down Tree
			strcat(info, "AHeader\n"); // Add "header"
			sprintf(work, "IOffset\t%08lXh\n", 0); strcat(info, work);
			sprintf(work, "ISize\t%ld\n", sizeof(DSF)); strcat(info, work);
#ifdef __GUNC__
			sprintf(work, "Ichunk_size\t%lld\n", dsf.chunk_size);
#else
			sprintf(work, "Ichunk_size\t%I64d\n", dsf.chunk_size);
#endif
			strcat(info, work);
#ifdef __GNUC__
			sprintf(work, "Ifile_size\t%lld\n", dsf.file_size);
#else
			sprintf(work, "Ifile_size\t%I64d\n", dsf.file_size);
#endif
			strcat(info, work);
#ifdef __GUNC__
			sprintf(work, "Iptr_meta\t%ld\n", dsf.ptr_meta);
#else
			sprintf(work, "Iptr_meta\t%I64d\n", dsf.ptr_meta);
#endif
			strcat(info, work);

			if (fseek_local(fp, dsf.chunk_size, SEEK_SET)) {
				break;
			}
			rs = fread(&dsf_fmt, 1, sizeof(DSF_FMT), fp);
			if (rs != sizeof(DSF_FMT)) {
				break;
			}
			strcat(info, "Afmt\n"); // Add "fmt"
			strcat(info, work);
			sprintf(work, "IOffset\t%08lXh\n", dsf.chunk_size);
			strcat(info, work);
			sprintf(work, "ISize\t%ld\n", sizeof(DSF_FMT));
			strcat(info, work);

			sprintf(work, "Ifmt_version\t%ld\n", dsf_fmt.fmt_version);
			strcat(info, work);
			sprintf(work, "Ifmt_id\t%ld\n", dsf_fmt.fmt_id);
			strcat(info, work);
			sprintf(work, "Ichannel_type\t%ld\n", dsf_fmt.channel_type);
			strcat(info, work);
			sprintf(work, "Ichannel_count\t%ld\n", dsf_fmt.channel_count);
			strcat(info, work);
			sprintf(work, "Isampling\t%ld\n", dsf_fmt.sampling);
			strcat(info, work);
			sprintf(work, "Isample_bit_count\t%ld\n",dsf_fmt.sample_bit_count); strcat(info, work);

#ifdef __GNUC__
			sprintf(work, "Isample_count\t%lld\n", dsf_fmt.sample_count);
#else
			sprintf(work, "Isample_count\t%I64d\n", dsf_fmt.sample_count);
#endif
			strcat(info, work);

			sprintf(work, "Iblock_size\t%ld\n", dsf_fmt.block_size);
			strcat(info, work);

			if (fseek_local(fp, dsf.chunk_size + dsf_fmt.chunk_size,SEEK_SET)) {
				break;
			}
			rs = fread(&dsf_data, 1, sizeof(DSF_DATA), fp);
			if (rs != sizeof(DSF_DATA)) {
				break;
			}

			strcat(info, "Adata\n"); // Add "data"
			strcat(info, work);
			sprintf(work, "IOffset\t%08lXh\n",dsf.chunk_size + dsf_fmt.chunk_size); strcat(info, work);
			sprintf(work, "ISize\t%ld\n", dsf_data.chunk_size + 12);
			strcat(info, work); retValue = STATUS_SUCCESS;
		}
		if (infoRF64File(filename, &inFmt, &inSample, NULL)) {
#if 0
			if (findWavChunk(fp, "fmt ", &chunkSize)) {
				strcat(info, "Afmt\n"); // Add "fmt"
				offset = ftell_local(fp);
				rs = fread(&wfex, 1, chunkSize, fp);
				if (rs != chunkSize) {
					break;
				}
				sprintf(work, "IID\tfmt \n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset);
				strcat(info, work); sprintf(work, "ISize\t%ld\n", rs);
				strcat(info, work);

				sprintf(work, "IwFormatTag\t%04x", wfex.Format.wFormatTag);
				strcat(info, work);
				if (wfex.Format.wFormatTag == WF_PCM) {
					strcat(info, "(WAVE_FORMAT_PCM)\n");
				} else if (wfex.Format.wFormatTag == WF_IEEE_FLOAT) {
					strcat(info, "(WAVE_FORMAT_IEEE_FLOAT)\n");
				} else if (wfex.Format.wFormatTag == 0xFFFE) {
					strcat(info, "(WAVE_FORMAT_EXTENSIBLE)\n");
				} else {
					strcat(info, "\n");
				}
				sprintf(work, "InChannels\t%d\n", wfex.Format.nChannels);
				strcat(info, work);
				sprintf(work, "InSamplesPerSec\t%ld\n",wfex.Format.nSamplesPerSec); strcat(info, work);
				sprintf(work, "IwBitsPerSample\t%d\n",wfex.Format.wBitsPerSample); strcat(info, work);
				sprintf(work, "InBlockAlign\t%d\n",wfex.Format.nBlockAlign); strcat(info, work);
				sprintf(work, "InAvgBytesPerSec\t%ld\n",wfex.Format.nAvgBytesPerSec); strcat(info, work);
				if (wfex.Format.wFormatTag == 0xFFFE) {
					sprintf(work, "IwValidBitsPerSample\t%d\n",wfex.Samples.wValidBitsPerSample);
					strcat(info, work);
					sprintf(work, "IwSamplesPerBlock\t%d\n",wfex.Samples.wSamplesPerBlock);
					strcat(info, work);
					sprintf(work, "IdwChannelMask\t%X\n",wfex.dwChannelMask);
					strcat(info, work);
				}
			}
			if (findWavChunk(fp, "data", &chunkSize)) {
				__int64 qw;
				strcat(info, "Adata\n"); // Add "data"
				offset = ftell_local(fp);
				sprintf(work, "IID\tdata\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset);
				strcat(info, work); qw = ds64.dataSizeHigh; qw <<= 32;
				qw |= ds64.dataSizeLow;
#ifdef __GNUC__
				sprintf(work, "ISize\t%lld\n", qw);
#else
				sprintf(work, "ISize\t%I64d\n", qw);
#endif
				strcat(info, work);

				if (qw > 0L) {
					qw /= (wfex.Format.wBitsPerSample / 8);
					qw /= (wfex.Format.nChannels);
#ifdef __GNUC__
					sprintf(work, "ISampleCount\t%lld\n", qw);
#else
					sprintf(work, "ISampleCount\t%I64d\n", qw);
#endif
					strcat(info, work);
				}
			}
			// List(info)
			if (findWavChunk(fp, "LIST", &chunkSize)) {
				if (chunkSize > 20) {
					offset = ftell_local(fp); listSize = chunkSize;
					listData = malloc(chunkSize);
					if (listData != NULL) {
						rs = fread(listData, 1, listSize, fp);
						if (rs != listSize) {
							free(listData);
							listData = NULL;
						} else {
							if (memcmp(listData, "INFO", 4) != 0) {
								free(listData);
								listData = NULL;
							}
						}
					}
				}
				if (listData != NULL) {
					strcat(info, "ALIST(INFO)\n");
					// Add "INFO"
					sprintf(work, "IID\tLIST(INFO)\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%ld\n", listSize + 8);
					strcat(info, work); p = listData + 4; listSize -= 4;
					while (listSize >= 8) {
						strcat(info, "I");
						memset(work, 0, 5);
						strncpy(work, p, 4);
						strcat(info, work);
						strcat(info, "\t");
						p += 4;
						len = *(long*)p;
						p += 4;
						memset(work, 0, 256);
						strncpy(work, p, len < 256 ? len : 255);
						if (len & 1) {
							len++;
						}
						strcat(info, work);
						strcat(info, "\n");
						p += len;
						if (listSize <= len + 8) {
							break;
						}
						listSize -= (len + 8);
					}
					free(listData);
					listData = NULL;
				}
			}

			if (findWavChunk(fp, "JUNK", &chunkSize)) {
				strcat(info, "AJUNK\n"); // Add "data"
				offset = ftell_local(fp);
				sprintf(work, "IID\tJUNK\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);
			}
			if (findWavChunk(fp, "ID3 ", &chunkSize)) {
				strcat(info, "AID3 \n"); // Add "data"
				offset = ftell_local(fp);
				sprintf(work, "IID\tID3 \n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);
			}
#endif
			retValue = STATUS_SUCCESS;
		}
	} while (0);

	if (listData != NULL) {
		free(listData);
	}
	if (bext != NULL) {
		free(bext);
	}
	if (fp != NULL) {
		fclose(fp);
	}

	return retValue; 
}
int findWavChunk(FILE * fp, char*id, long*size)
{
	WAVE_HEADER *wh;
	long rs;
	__int64 dataSize;
	__int64 seekPtr;
	int flg_datachunk;
	DWORD dataLow,dataHigh;
	int found;
	wh = (WAVE_HEADER*)tempBuffer;

	found = 0;
	flg_datachunk = 0;
	do {
		// WAV header
		if (fseek_local(fp, 0, SEEK_SET)) {
			break;
		}
		rs = fread(tempBuffer, 1, 12, fp);
		if (rs != 12) {
			break;
		}
		if (!(memcmp(tempBuffer, "RIFF", 4) == 0 || memcmp(tempBuffer,"riff", 4) == 0 || 
			memcmp(tempBuffer, "RF64",4) == 0 || memcmp(tempBuffer, "rf64", 4) == 0)) {
			break;
		}
		if (!(memcmp(tempBuffer + 8, "WAVE",4) == 0 || memcmp(tempBuffer + 8, "wave", 4) == 0)) {
			break;
		}
		seekPtr = 12;
		// 指定チャンクのサーチ
		for (; ; ) {
			if (fseek_local(fp, seekPtr, SEEK_SET)) {
				break;
			}
			rs = fread(wh, 8, 1, fp);
			if (rs != 1) {
				break;
			}
			if (wh->id[0] == 'd' && wh->id[1] == 's' && wh->id[2] == '6' && wh->id[3] == '4') {
				if (fread(&dataLow,1,4,fp) != 4) {
					break;
				}
				if (fread(&dataHigh,1,4,fp) != 4) {
					break;
				}
				if (fread(&dataLow,1,4,fp) != 4) {
					break;
				}
				if (fread(&dataHigh,1,4,fp) != 4) {
					break;
				}
				if (fseek_local(fp, seekPtr + 8, SEEK_SET)) {
					break;
				}
				dataSize = dataLow;dataSize <<= 32;dataSize |= dataHigh;
			}
			if (wh->id[0] == 'd' && wh->id[1] == 'a' && wh->id[2] == 't' && wh->id[3] == 'a') {
				flg_datachunk = 1;
			}
			if (memcmp(wh->id, id, 4) == 0 && wh->size > 0) {
				found = 1;
				*size = wh->size;
				break;
			}
			if (flg_datachunk == 0) {
				if (wh->size & 1) {
					// 奇数なら00のパディングがある
					seekPtr += wh->size + 9;
				} else {
					seekPtr += wh->size + 8;
				}
			} else {
				flg_datachunk = 0;
				if (wh->size != 0xFFFFFFFF) {
					dataSize = wh->size;
				}
				seekPtr += dataSize + 8;
			}
		}
	}
	while (0);
	return found;
}
int findWavChunk_N(FILE * fp, int n,char *tag,long*size)
{
	WAVE_HEADER *wh;
	long rs;
	__int64 dataSize;
	__int64 seekPtr;
	int flg_datachunk;
	DWORD dataLow,dataHigh;
	int found;
	int i;
	wh = (WAVE_HEADER*)tempBuffer;

	found = 0;
	flg_datachunk = 0;
	do {
		// WAV header
		if (fseek_local(fp, 0, SEEK_SET)) {
			break;
		}
		rs = fread(tempBuffer, 1, 12, fp);
		if (rs != 12) {
			break;
		}
		if (!(memcmp(tempBuffer, "RIFF", 4) == 0 || memcmp(tempBuffer,"riff", 4) == 0 || 
			memcmp(tempBuffer, "RF64",4) == 0 || memcmp(tempBuffer, "rf64", 4) == 0)) {
			break;
		}
		if (!(memcmp(tempBuffer + 8, "WAVE",4) == 0 || memcmp(tempBuffer + 8, "wave", 4) == 0)) {
			break;
		}
		seekPtr = 12;
		// 指定チャンクのサーチ
		for (i = 0;i < n;) {
			i++;
			if (fseek_local(fp, seekPtr, SEEK_SET)) {
				break;
			}
			rs = fread(wh, 8, 1, fp);
			if (rs != 1) {
				break;
			}
if (1) {
	char ssss[256];
	sprintf(ssss,"%c%c%c%c:%08lx",wh->id[0],wh->id[1],wh->id[2],wh->id[3],seekPtr);
	PRINT_LOG(ssss);
}
			if (wh->id[0] == 'd' && wh->id[1] == 's' && wh->id[2] == '6' && wh->id[3] == '4') {
				if (fread(&dataLow,1,4,fp) != 4) {
					break;
				}
				if (fread(&dataHigh,1,4,fp) != 4) {
					break;
				}
				if (fread(&dataLow,1,4,fp) != 4) {
					break;
				}
				if (fread(&dataHigh,1,4,fp) != 4) {
					break;
				}
				if (fseek_local(fp, seekPtr + 8, SEEK_SET)) {
					break;
				}
				dataSize = dataHigh;dataSize <<= 32;dataSize |= dataLow;
			}
			if (wh->id[0] == 'd' && wh->id[1] == 'a' && wh->id[2] == 't' && wh->id[3] == 'a') {
				flg_datachunk = 1;
			}
			if (i == n && wh->size > 0) {
				found = 1;
				tag[0] = wh->id[0];
				tag[1] = wh->id[1];
				tag[2] = wh->id[2];
				tag[3] = wh->id[3];
				*size = wh->size;
if (1) {
	char ssss[256];
	sprintf(ssss,"Found:%c%c%c%c:%08lx",wh->id[0],wh->id[1],wh->id[2],wh->id[3],seekPtr);
	PRINT_LOG(ssss);
}
				break;
			}
			if (flg_datachunk == 0) {
				if (wh->size & 1) {
					// 奇数なら00のパディングがある
					seekPtr += wh->size + 9;
				} else {
					seekPtr += wh->size + 8;
				}
			} else {
				flg_datachunk = 0;
				if (wh->size != 0xFFFFFFFF) {
					dataSize = wh->size;
				}
				seekPtr += dataSize + 8;
			}
		}
	}
	while (0);
	return found;
}

static __int64 ftell_local(FILE * fp) {
	__int64 pos;
#ifdef __GNUC__
	pos = ftello64(fp);
#else
	pos = _ftelli64(fp);
#endif
	return pos;
}
static int fseek_local(FILE * fp, __int64 offset,int origin) {
	int ret;
#ifdef __GNUC__
	ret = fseeko64(fp, offset, origin);
#else
	ret = _fseeki64(fp, offset, origin);
#endif
	return ret;
}
#ifndef _Windows
char *utf8_to_utf16(char *u8_str,int *to_size)
{
	char *to_str;
	size_t u8_size;
	size_t u16_size;
	char *p_in,*p_out;
	size_t size_in,size_out;
	iconv_t ic;
	int retry;
	int mc;
	int j;

	u8_size = strlen(u8_str);
	to_str = NULL;

	ic = iconv_open("UTF-16LE","UTF-8");
	if (ic != (iconv_t)-1) {
		for (retry = 5,mc = 3;retry >= 0;retry--,mc++) {
			if (to_str != NULL) {free(to_str);to_str = NULL;}
			to_str = malloc((u8_size * mc) + 4);
			if (to_str != NULL) {
				size_t r;
				p_in  = u8_str;
				p_out = to_str + 2;
				size_in  = u8_size;
				size_out = (u8_size * mc);
				errno = 0;
				r = iconv(ic,&p_in,&size_in,&p_out,&size_out);
				if (r != (iconv_t)(-1)) {
					// 正常終了
					iconv_close(ic);
					p_out[0] = '\0';
					p_out[1] = '\0';
					if (to_str[2] == 0xff && to_str[3] == 0xfe) {
						*to_size = u16_strlen(to_str + 2);
						for (j = 0;j < *to_size;j++) {
							to_str[j] = to_str[j + 2];
						}
					} else {
						to_str[0] = 0xff;
						to_str[1] = 0xfe;
						*to_size = u16_strlen(to_str);
					}
					return to_str;
				} else {
//					fprintf(stdout,"iconv error:%d,%s\n",errno,strerror(errno));
				}
			}
		}
		iconv_close(ic);
	}
	if (to_str != NULL) {free(to_str);to_str = NULL;}
	to_str = malloc(u8_size + 1);
	if (to_str != NULL) {
		strcpy(to_str,u8_str);
		*to_size = u8_size + 1;
	}
	return to_str;
}

char *utf8_to_sjis(char *u8_str,int *to_size)
{
	char *to_str;
	size_t u8_size;
	size_t sjis_size;
	char *p_in,*p_out;
	size_t size_in,size_out;
	iconv_t ic;
	int retry;
	int mc;

	u8_size = strlen(u8_str);
	to_str = NULL;

	ic = iconv_open("CP932","UTF-8");
	if (ic != (iconv_t)-1) {
		for (retry = 5,mc = 3;retry >= 0;retry--,mc++) {
			if (to_str != NULL) {free(to_str);to_str = NULL;}
			to_str = malloc((u8_size * mc) + 1);
			if (to_str != NULL) {
				size_t r;
				p_in  = u8_str;
				p_out = to_str;
				size_in  = u8_size;
				size_out = (u8_size * mc);
				errno = 0;
				r = iconv(ic,&p_in,&size_in,&p_out,&size_out);
				if (r != (iconv_t)(-1)) {
					// 正常終了
					iconv_close(ic);
					p_out[0] = '\0';
					*to_size = strlen(to_str) + 1;
					return to_str;
				} else {
//					fprintf(stdout,"iconv error:%d,%s\n",errno,strerror(errno));
				}
			}
		}
		iconv_close(ic);
	}
	if (to_str != NULL) {free(to_str);to_str = NULL;}
	to_str = malloc(u8_size + 1);
	if (to_str != NULL) {
		strcpy(to_str,u8_str);
		*to_size = u8_size + 1;
	}
	return to_str;
}
static int u16_strlen(char *u16_str)
{
	int count;
	WORD *wp = (WORD*)u16_str;
	for (count = 0;*wp != 0;wp++,count++)
		;

	return (count + 1) * 2;
}
#else
char *utf8_to_utf16(char *u8_str,int *to_size)
{
}
char *utf8_to_sjis(char *u8_str,int *to_size)
{
}
#endif
