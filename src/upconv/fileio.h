/****************************************************************************/
/* fileio (C) 2011 By 59414d41												*/
/* テンポラリファイル入出力IO												*/
/*																			*/
/****************************************************************************/

#ifndef FILEIO_H
#define FILEIO_H
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>

typedef signed long long fio_size;

#define FIO_MODE_R		(1)				// 読み取り専用
#define FIO_MODE_W		(2)				// 書き込み専用

#define FIO_ERR_READ	(-1)
#define FIO_ERR_WRITE	(-2)
#define FIO_ERR_MEMORY	(-3)
#define FIO_ERR_OTHER	(-4)

typedef struct {
	FILE *fp;
	fio_size data_size;
	fio_size data_offset;
	fio_size data_maxsize;
	void *p_data;
	int magic;
	int	mode;
	int error;
	int  block_count;
	int	 block_max;
	long block_size;
	unsigned char *ibuffer;
	int ibuffer_size;
	char fname[_MAX_PATH];
} FIO;

void fio_open(FIO *fio,char *filename,int mode);
void fio_close(FIO *fio);
void fio_get_filesize(FIO *fio,fio_size *size);
void fio_get_datasize(FIO *fio,fio_size *size);
void fio_set_tmpfile(FIO *fio,char *filename);
void fio_set_memory_limit(FIO *fio,int block_sz,int max);
void fio_set_maxsize(FIO *fio,fio_size max);
void fio_setmode_r(FIO *fiow,FIO *fior,char *filename);
void fio_seek(FIO *fio,fio_size offset,int orign);
fio_size fio_read(void *buf,fio_size size,fio_size n,FIO *fio);
fio_size fio_write(void *buf,fio_size size,fio_size n,FIO *fio);
void fio_flush(FIO *fio);
void fio_rewind(FIO *fio);
fio_size fio_tell(FIO *fio);
#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif

