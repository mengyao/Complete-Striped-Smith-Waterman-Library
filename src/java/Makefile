CC = gcc
CFLAGS := -Wall -O0 -g -pipe #-pg
JNI_JAR = ssw.jar
JNI_LIB = libsswjni.so
INLCUDES=  -I../ -I"$(JAVA_HOME)/include" -I"$(JAVA_HOME)/include/linux"
JNI_JAVA_OBJ = ssw/Aligner.class ssw/Alignment.class ssw/Example.class

.PHONY: all clean

all: $(JNI_JAR) $(JNI_LIB) 

$(JNI_LIB): sswjni.c ../ssw.c ../ssw.h
	$(CC) $(CFLAGS) $(INLCUDES) -fPIC -shared -rdynamic -o $@ $< ../ssw.c 

$(JNI_JAR): $(JNI_JAVA_OBJ)
	jar cvfe $@ ssw.Example $^

%.class: %.java
	javac $<

clean:
	-rm -f $(JNI_LIB) $(JNI_JAR) $(JNI_JAVA_OBJ)
