#include <genmap-impl.h>

#define HEADER_LEN 132

int parRSBWriteMapFile(int nelt,int nv,int *pmap,int *vtx,char *fileName,
		       MPI_Comm comm){
  char *version="#v001";
  int nelgt=nelt;
  MPI_Allreduce(&nelt,&nelgt,1,MPI_INT,MPI_SUM,comm);
  const int npts=nelgt*nv;
  const int depth=(int)log2(1.0*nelgt);
  const int d2=(int)(pow(2,depth)+0.5);
  int nactive=nelgt,nrnk=nelgt,noutflow=0;

  char header[HEADER_LEN];
  header[HEADER_LEN]='\0';
  sprintf(header,"%5s%12d%12d%12d%12d%12d%12d%12d",version,
		  nelgt,nactive,depth,d2,npts,nrnk,noutflow);
  memset(header+strlen(header),' ',HEADER_LEN-strlen(header));
  header[HEADER_LEN]='\0';

  MPI_Info infoIn;
  MPI_Info_create(&infoIn);
  MPI_Info_set(infoIn,"access_style","write_once,random");

  int errs=0;
  MPI_File file;
  int err=MPI_File_open(comm,fileName,
		MPI_MODE_WRONLY|MPI_MODE_CREATE,
		infoIn,&file);
  if(err){
    errs++;
    MPI_Abort(comm,911);
  }

  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  printf("nelt: %d rank=%d\n",nelt,rank);
  int writeSize=0;
  if(rank==0)
    writeSize=HEADER_LEN*sizeof(char)+sizeof(float);
  writeSize+=(nv+1)*nelt*sizeof(int);

  char *buf=(char*)malloc(writeSize),*buf0=buf;

  float test=6.54321;
  if(rank==0){
    memcpy(buf0,header,HEADER_LEN*sizeof(char)),buf0+=HEADER_LEN;
    memcpy(buf0,&test,sizeof(float)),buf0+=sizeof(float);
  }

  int i;
  for(i=0;i<nelt;i++){
    memcpy(buf0,&pmap[i],sizeof(int)),buf0+=sizeof(int);
    memcpy(buf0,&vtx[i*nv],sizeof(int)*nv),buf0+=nv*sizeof(int);
  }

  MPI_Status status;
  err=MPI_File_write_ordered(file,buf,writeSize,MPI_BYTE,&status);
  if(err) errs++;

#if 0
  MPI_Info_set(infoIn,"access_style","read_once");
  MPI_File_seek_shared(file,0,MPI_SEEK_SET);

  MPI_File_set_info(file,infoIn);
  MPI_Info_free(&infoIn);
  buf[0]=-1;

  err=MPI_File_read_ordered(file,buf,writeSize,MPI_BYTE,&status);
  if(err) errs++;
  int count;
  MPI_Get_count(&status,MPI_BYTE,&count);
  if(count!=writeSize){
    errs++;
    printf("Expected to read one int, read %d\n ints",count);
  }
#endif

  err=MPI_File_close(&file);
  if(err) errs++;

  MPI_Barrier(comm);
  free(buf);

  return errs;
}
