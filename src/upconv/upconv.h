#ifndef UPCONV_H
#define UPCONV_H

#ifndef BYTE
typedef unsigned char 	BYTE;
#endif
#ifndef WORD
typedef unsigned short	WORD;
#endif
#ifndef DWORD
typedef unsigned long	DWORD;
#endif

/* Return Status */
#define STATUS_SUCCESS			(0)		/* 正常終了 */
#define STATUS_NOT_IMPLEMENT	(-1)	/* その機能はインプリメントされていない */
#define STATUS_UNKNOWN_FORMAT	(-2)	/* 未知のフォーマット */
#define STATUS_DATA_ERR			(-3)	/* データが壊れている */
#define STATUS_MEM_ALLOC_ERR	(-4)	/* メモリーが確保出来ない */
#define STATUS_PARAMETER_ERR	(-5)
#define STATUS_FILE_READ_ERR	(-6)	/* ファイルリードエラー */
#define STATUS_FILE_WRITE_ERR	(-7)	/* ファイルライトエラー */
#define STATUS_EXEC_FAIL		(-8)	/* プロセス実行エラー */
#define STATUS_PLUGIN_ERR		(-15)	/* 内部エラー */

typedef struct {
	int		state;
	HANDLE	ps;
	HANDLE	thread;
	HANDLE	hStdOutRead;
} STARTEXEC_INFO;

#endif

