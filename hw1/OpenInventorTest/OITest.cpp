/*
 * OITest.c
 *
 *  Created on: Sep 15, 2018
 *      Author: joe
 */
#include <Inventor/SoOutput.h>
#include <Inventor/SoDB.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/nodes/SoCone.h>
#include <Inventor/nodes/SoSeparator.h>
#include <stdlib.h>
#include <stdio.h>

static char * buffer;
static size_t buffer_size = 0;

static void *
buffer_realloc(void * bufptr, size_t size)
{
  buffer = (char *)realloc(bufptr, size);
  buffer_size = size;
  return buffer;
}

static SbString
buffer_writeaction(SoNode * root)
{
  SoOutput out;
  buffer = (char *)malloc(1024);
  buffer_size = 1024;
  out.setBuffer(buffer, buffer_size, buffer_realloc);

  SoWriteAction wa(&out);
  wa.apply(root);

  SbString s(buffer);
  free(buffer);
  return s;
}

int
main(int argc, char ** argv)
{
  SoDB::init();

  SoSeparator * root = new SoSeparator;
  root->ref();

  root->addChild(new SoCone);

  SbString s = buffer_writeaction(root);
  (void)fprintf(stdout, "%s\n", s.getString());

  root->unref();
  return 0;
}
