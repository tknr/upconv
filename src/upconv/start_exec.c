#include <stdio.h>
#include <string.h>
#include <windows.h>

#include "upconv.h"

int start_exec(int argc,char *argv[],int cpu_pri,HANDLE *ret_hStdOutRead,HANDLE *ret_process_id,HANDLE *ret_thread_id)
{
	HANDLE hStdOutReadTmp,hStdOutRead,hStdOutWrite;
	HANDLE hStdErrorWrite;
	SECURITY_ATTRIBUTES sa;
	PROCESS_INFORMATION pi;
	STARTUPINFO si;
	int i;
	char cmd[4096];
	char work[2048];

	memset(&sa,0,sizeof(SECURITY_ATTRIBUTES));
	sa.nLength= sizeof(SECURITY_ATTRIBUTES);
	sa.bInheritHandle = TRUE;

	if (ret_hStdOutRead != NULL) {
		// Read,Write用のパイプを作成する
		if (!CreatePipe(&hStdOutReadTmp,&hStdOutWrite,&sa,0)) return -1;

		// 標準エラー出力を複製する
		if (!DuplicateHandle(GetCurrentProcess(),hStdOutWrite,GetCurrentProcess(),&hStdErrorWrite,0,TRUE,DUPLICATE_SAME_ACCESS)) return -1;

		if (!DuplicateHandle(GetCurrentProcess(),hStdOutReadTmp,GetCurrentProcess(),&hStdOutRead,0,FALSE,DUPLICATE_SAME_ACCESS)) return -1;

		CloseHandle(hStdOutReadTmp);
	}

	memset(&si,0,sizeof(STARTUPINFO));
	si.cb = sizeof(STARTUPINFO);
	si.dwFlags = STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW | 0x00000020;
	si.wShowWindow = SW_HIDE;
	si.hStdInput  = GetStdHandle(STD_INPUT_HANDLE);
	if (ret_hStdOutRead != NULL) {
		si.hStdOutput = hStdOutWrite;
		si.hStdError  = hStdErrorWrite;
	} else {
		si.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
		si.hStdError  = GetStdHandle(STD_ERROR_HANDLE);
	}

	// コマンドの組み立て
	cmd[0] = '\0';
	for (i = 0;i < argc;i++) {
		if (argv[i][0] != '|') {
			sprintf(work,"\"%s\" ",argv[i]);
		} else {
			sprintf(work,"%s ",&argv[i][1]);
		}
		strcat(cmd,work);
	}

	if (!CreateProcess(NULL,cmd,NULL,NULL,TRUE,/* CREATE_NEW_CONSOLE | */ IDLE_PRIORITY_CLASS /* | 0x00000200 */,NULL,NULL,&si,&pi)) return -1;

	if (ret_hStdOutRead != NULL) *ret_hStdOutRead = hStdOutRead;
	if (ret_process_id  != NULL) *ret_process_id = pi.hProcess;
	if (ret_thread_id   != NULL) *ret_thread_id   = pi.hThread;
	if (cpu_pri == 0) {
		SetPriorityClass(pi.hProcess,IDLE_PRIORITY_CLASS);
	} else {
		SetPriorityClass(pi.hProcess,NORMAL_PRIORITY_CLASS);
	}
	if (ret_hStdOutRead != NULL) {
		CloseHandle(hStdOutWrite);
		CloseHandle(hStdErrorWrite);
	}
	
	return 0;
}
int finish_exec(HANDLE *ret_hStdOutRead,HANDLE *ret_process_id,HANDLE *ret_thread_id)
{
	if (ret_hStdOutRead != NULL && *ret_hStdOutRead != INVALID_HANDLE_VALUE) CloseHandle(*ret_hStdOutRead);
	if (ret_process_id  != NULL && *ret_process_id  != INVALID_HANDLE_VALUE) CloseHandle(*ret_process_id);
	if (ret_thread_id  != NULL && *ret_thread_id  != INVALID_HANDLE_VALUE) CloseHandle(*ret_thread_id);

	return 0;
}

