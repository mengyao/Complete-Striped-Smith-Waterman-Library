#include <jni.h>
#include <stdio.h>
#include "../ssw.h"

jstring s_alignToJniCigarString(JNIEnv* env, s_align* align) {
	if (align->cigar != NULL && align->cigarLen > 0) {
		char* buffer = (char *)malloc(16 * align->cigarLen + 1);
		buffer[0] = '\0';
		char* currentBufferPosition = buffer;
		int i;
		for (i = 0; i < align->cigarLen; i++) {
			int charsPrinted = sprintf(currentBufferPosition, "%d%c", cigar_int_to_len(align->cigar[i]), cigar_int_to_op(align->cigar[i]));
			currentBufferPosition += charsPrinted;
			printf("Cigar is now %s\n", buffer);
		}
		jstring jstrBuf = (*env)->NewStringUTF(env, buffer);
		free(buffer);
		return jstrBuf;
	}
	return NULL;
}
int maxArrayValue(jbyte* array, int len) {
	int i, max = 0;
	for (i = 0; i < len; i++) {
		if (array[i] > max) {
			max = array[i];
		}
	}
	return max;
}
void debugPrintjbyteArray(jbyte* ptr, jsize len) {
	int i;
	printf("{ ");
	for (i = 0; i < len; i++) {
		printf("%d,", ptr[i]);
	}
	printf("} ");
}

JNIEXPORT jlong JNICALL Java_ssw_Aligner_initprofile(JNIEnv *env, jclass cls, jbyteArray read, jbyteArray matrix, jint matrixSize, jint score_size) {
	s_profile* profile = NULL;
	jbyte* readPtr = (*env)->GetByteArrayElements(env, read, NULL);
	jsize readLen = (*env)->GetArrayLength(env, read);
	jbyte* matrixPtr = (*env)->GetByteArrayElements(env, matrix, NULL);
	jsize matrixLen = (*env)->GetArrayLength(env, matrix);
	printf("Java_ssw_Aligner_initprofile: "); debugPrintjbyteArray(readPtr, readLen) ; debugPrintjbyteArray(matrixPtr, matrixLen) ; printf(", %d, %d\n", matrixSize, score_size); fflush(stdout);
	if (matrixLen == matrixSize * matrixSize && maxArrayValue(readPtr, readLen) < matrixSize) {
		profile = ssw_init(readPtr, readLen, matrixPtr, matrixSize, (int8_t)score_size);
	}
	(*env)->ReleaseByteArrayElements(env, read, readPtr, JNI_ABORT);
	(*env)->ReleaseByteArrayElements(env, matrix, matrixPtr, JNI_ABORT);
	return (jlong)profile;
}
JNIEXPORT void JNICALL Java_ssw_Aligner_destroyprofile(JNIEnv *env, jclass cls, jlong profilePtr) {
	init_destroy((s_profile*)profilePtr);
}
JNIEXPORT jobject JNICALL Java_ssw_Aligner_align(JNIEnv *env, jclass cls, jlong profilePtr,
		jbyteArray ref,
		jint gapOpen,
		jint gapExtend,
		jint flag,
		jshort filters,
		jint filterd,
		jint maskLen
		) {
	jobject result = NULL;
	s_align* align = NULL;
	jbyte* refPtr = (*env)->GetByteArrayElements(env, ref, NULL);
	jsize refLen = (*env)->GetArrayLength(env, ref);
	printf("Java_ssw_Aligner_align: %p", (void*)profilePtr); debugPrintjbyteArray(refPtr, refLen) ; printf(", %d, %d, 0x%02x, %d, %d, %d, \n", gapOpen, gapExtend, flag, filters, filterd, maskLen); fflush(stdout);
	align = ssw_align((s_profile*)profilePtr, refPtr, refLen, gapOpen, gapExtend, flag, filters, filterd, maskLen);
	(*env)->ReleaseByteArrayElements(env, ref, refPtr, JNI_ABORT);
	if (align != NULL) {
		jclass clazz = (*env)->FindClass(env, "ssw/Alignment");
		jmethodID constructor = (*env)->GetMethodID(env, clazz, "<init>", "(SSIIIIILjava/lang/String;)V");
		result = (*env)->NewObject(env, cls, constructor,
			align->score1,
			align->score2,
			align->ref_begin1,
			align->ref_end1,
			align->read_begin1,
			align->read_end1,
			align->ref_end2,
			s_alignToJniCigarString(env, align));
		align_destroy(align);
	}
	return result;
}


