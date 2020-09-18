!To compile, run the following line (without leading exclamation mark):
!gfortran ASoT_Functions_Fortran.F90 SurfaceLine_Advance.F90 -o SurfaceLine_Advance

!Filenames, directories, lines in control files are capped at 200 characters by len=200 magic number

program SurfaceLine_Advance
      
      use ASoT_Functions_Fortran
      
      implicit none
      
      real(4),parameter :: PI=4.0*atan(1.0)
      character (len = 200) :: CNT_file      !Control file for this program
      character (len = 200) :: ARMS_CNT_file !For the ARMS control file
      character (len = 200) :: SurfaceLine_initial
      character (len = 200) :: SurfaceLine_current
      real(4) :: t_start
      real(4) :: delta_t
      real(4) :: solar_Radius
      integer(4) :: num_t
      real(4),allocatable :: plot_t(:)
      character,allocatable :: SurfaceLine_outputs(:,:)
      integer(4) :: num_static_cage
      integer(4) :: num_dynamic_cage
      integer(4),allocatable :: len_static_cages(:)
      integer(4),allocatable :: len_dynamic_cages(:)
      real(4),allocatable :: static_cages(:,:)
      real(4),allocatable :: dynamic_cages(:,:)
      real(4),allocatable :: cage_out(:,:,:)
      
      integer(4) :: idx_c,idx_f,idx,idx2,sum_static_cage,static_cage_st,in_cage_st
      
      CNT_file='SurfaceLine_Advance.cnt'
      call parse_SurfaceLine_CNT(CNT_file,ARMS_CNT_file,SurfaceLine_initial,t_start,delta_t,solar_Radius,&
      &num_t,plot_t,SurfaceLine_outputs)
      
      call surfaceline_read_binary_cage(SurfaceLine_initial,num_static_cage,len_static_cages,num_dynamic_cage,&
      &len_dynamic_cages,static_cages,dynamic_cages)
      
      do idx_f=1,num_t
	do idx=1,200
	  SurfaceLine_current(idx:idx)=SurfaceLine_outputs(idx,idx_f)
	end do
        sum_static_cage=0
        open (unit=6,file=SurfaceLine_current,form='unformatted',status='REPLACE',recl=4)
	close(6)
        open (unit=6,file=SurfaceLine_current,form='unformatted',access='direct',recl=4)
        write (6,rec=1) num_static_cage
	do idx=1,num_static_cage
	  write (6,rec=idx+1) len_static_cages(idx)
	  sum_static_cage=sum_static_cage+len_static_cages(idx)
	end do
	
        write (6,rec=2+num_static_cage) num_dynamic_cage
	do idx=1,num_dynamic_cage
	  write (6,rec=idx+2+num_static_cage) len_dynamic_cages(idx)
	end do
	
	static_cage_st=0
	do idx=1,num_static_cage
	  do idx2=1,len_static_cages(idx)
	    write (6,rec=idx2+static_cage_st*2+2+num_static_cage+num_dynamic_cage) static_cages(1,idx2+static_cage_st)
	  end do
	  do idx2=1,len_static_cages(idx)
	    write (6,rec=idx2+static_cage_st*2+2+num_static_cage+num_dynamic_cage+len_static_cages(idx))&
	    & static_cages(2,idx2+static_cage_st)
	  end do
	  static_cage_st=static_cage_st+len_static_cages(idx)
	end do
	close(6)
      end do
      
      in_cage_st=0
      do idx_c=1,num_dynamic_cage
        call surfaceline_advance_cage(len_dynamic_cages(idx_c),dynamic_cages(:,in_cage_st+1:&
        &in_cage_st+len_dynamic_cages(idx_c)),t_start,num_t,plot_t,delta_t,ARMS_CNT_file,solar_Radius,cage_out)
        
        do idx_f=1,num_t
	  do idx=1,200
	    SurfaceLine_current(idx:idx)=SurfaceLine_outputs(idx,idx_f)
	  end do
          open (unit=6,file=SurfaceLine_current,form='unformatted',access='direct',recl=4)

	  do idx2=1,len_dynamic_cages(idx_c)
	    write (6,rec=idx2+in_cage_st*2+2+num_static_cage+num_dynamic_cage+2*sum_static_cage)&
	    & cage_out(1,idx2,idx_f)
	  end do
	  do idx2=1,len_dynamic_cages(idx_c)
	    write (6,rec=idx2+in_cage_st*2+2+num_static_cage+num_dynamic_cage+2*sum_static_cage&
	    &+len_dynamic_cages(idx_c)) cage_out(2,idx2,idx_f)
	  end do
	  close(6)
        end do
        in_cage_st=in_cage_st+len_dynamic_cages(idx_c)

        
        deallocate(cage_out)
      end do
      
end program SurfaceLine_Advance


