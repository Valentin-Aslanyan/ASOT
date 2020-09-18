!Filenames, directories, lines in control files are capped at 200 characters by len=200 magic number

module ASoT_Functions_Fortran
      implicit none

      public :: ARMS_Gauss,ARMS_DGauss,surface_move_rk_4,surfaceline_read_binary_cage
      public :: surfaceline_write_binary_cage,parse_SurfaceLine_CNT,str_raise,str_lower
      public :: string_in_string,string_stats_letter,string_split_num,surfaceline_advance_cage
      public :: parse_QSL_Rbinfile

      contains


      function ARMS_Gauss(x,x_c,c) result(res)
        real(4), intent(in) :: x
        real(4), intent(in) :: x_c
        real(4), intent(in) :: c
        real(4) :: res
        res=EXP(-c*(x-x_c)*(x-x_c))
      end function ARMS_Gauss


      function ARMS_DGauss(x,x_c,c) result(res)
        real(4), intent(in) :: x
        real(4), intent(in) :: x_c
        real(4), intent(in) :: c
        real(4) :: res
        res=-2.0*c*(x-x_c)*EXP(-c*(x-x_c)*(x-x_c))
      end function ARMS_DGauss


      function surface_move_rk_4(R,theta,phi,t,delta_t,num_prof,Time_profiles,sflow_parameters) result(arr_out)
        real(4), intent(in) :: R
        real(4), intent(in) :: theta
        real(4), intent(in) :: phi
        real(4), intent(in) :: t
        real(4), intent(in) :: delta_t
        integer(4), intent(in) :: num_prof
        integer(4), intent(in) :: Time_profiles(:)
        real(4), intent(in) :: sflow_parameters(:,:)
        real(4) :: arr_out(2)
        
        integer(4) :: idx_prof
        real(4) :: k_1_th,k_1_ph,k_2_th,k_2_ph,k_3_th,k_3_ph,k_4_th,k_4_ph
        real(4) :: A_time,t_p,v_theta,v_phi
        real(4) :: sfl2l,sfl2r,sfl2c,k2sfl,v0sfl2,sfl3l,sfl3r,sfl3c,k3sfl,v0sfl3,tlsfl,trsfl,tcsfl,ktsfl
        real(4),parameter :: PI=4.0*atan(1.0)
        
        
        k_1_th=0.0
        k_1_ph=0.0
        k_2_th=0.0
        k_2_ph=0.0
        k_3_th=0.0
        k_3_ph=0.0
        k_4_th=0.0
        k_4_ph=0.0
        do idx_prof=1,num_prof
          if (Time_profiles(idx_prof)>1) then
            sfl2l=sflow_parameters(idx_prof,6)
            sfl2r=sflow_parameters(idx_prof,7)
            sfl2c=sflow_parameters(idx_prof,8)
            k2sfl=sflow_parameters(idx_prof,9)
            v0sfl2=sflow_parameters(idx_prof,10)
            sfl3l=sflow_parameters(idx_prof,11)
            sfl3r=sflow_parameters(idx_prof,12)
            sfl3c=sflow_parameters(idx_prof,13)
            k3sfl=sflow_parameters(idx_prof,14)
            v0sfl3=sflow_parameters(idx_prof,15)
            tlsfl=sflow_parameters(idx_prof,16)
            trsfl=sflow_parameters(idx_prof,17)
            tcsfl=sflow_parameters(idx_prof,18)
            ktsfl=sflow_parameters(idx_prof,19)
            if (t>=tlsfl .and. t<=trsfl) then
              A_time=0.0
              if (Time_profiles(idx_prof)==2) then
                t_p=2.0*PI*ktsfl/(trsfl-tlsfl)
                A_time=0.5*(1.0-cos(t_p*(t-tcsfl)))
              else if (Time_profiles(idx_prof)==3) then
                t_p=2.0*PI*ktsfl/(trsfl-tlsfl)
                A_time=sin(t_p*(t-tcsfl))
              end if
              v_theta=v0sfl2*ARMS_Gauss(theta,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))&
              &*ARMS_DGauss(phi,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/sin(PI*theta)/PI*A_time
              v_phi=v0sfl3*ARMS_DGauss(theta,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))&
              &*ARMS_Gauss(phi,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/PI*A_time
              k_1_th=k_1_th+v_theta/R/PI
              k_1_ph=k_1_ph+v_phi/R/PI/sin(PI*theta)
            end if
          end if
        end do
        do idx_prof=1,num_prof
          if (Time_profiles(idx_prof)>1) then
            sfl2l=sflow_parameters(idx_prof,6)
            sfl2r=sflow_parameters(idx_prof,7)
            sfl2c=sflow_parameters(idx_prof,8)
            k2sfl=sflow_parameters(idx_prof,9)
            v0sfl2=sflow_parameters(idx_prof,10)
            sfl3l=sflow_parameters(idx_prof,11)
            sfl3r=sflow_parameters(idx_prof,12)
            sfl3c=sflow_parameters(idx_prof,13)
            k3sfl=sflow_parameters(idx_prof,14)
            v0sfl3=sflow_parameters(idx_prof,15)
            tlsfl=sflow_parameters(idx_prof,16)
            trsfl=sflow_parameters(idx_prof,17)
            tcsfl=sflow_parameters(idx_prof,18)
            ktsfl=sflow_parameters(idx_prof,19)
            if (t+0.5*delta_t>=tlsfl .and. t+0.5*delta_t<=trsfl) then
              A_time=0.0
              if (Time_profiles(idx_prof)==2) then
                t_p=2.0*PI*ktsfl/(trsfl-tlsfl)
                A_time=0.5*(1.0-cos(t_p*(t+0.5*delta_t-tcsfl)))
              else if (Time_profiles(idx_prof)==3) then
                t_p=2.0*PI*ktsfl/(trsfl-tlsfl)
                A_time=sin(t_p*(t+0.5*delta_t-tcsfl))
              end if
              v_theta=v0sfl2*ARMS_Gauss(theta+0.5*delta_t*k_1_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))&
              &*ARMS_DGauss(phi+0.5*delta_t*k_1_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))&
              &/sin(PI*(theta+0.5*delta_t*k_1_th))/PI*A_time
              v_phi=v0sfl3*ARMS_DGauss(theta+0.5*delta_t*k_1_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))&
              &*ARMS_Gauss(phi+0.5*delta_t*k_1_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/PI*A_time
              k_2_th=k_2_th+v_theta/R/PI
              k_2_ph=k_2_ph+v_phi/R/PI/sin(PI*(theta+0.5*delta_t*k_1_th))
            end if
          end if
        end do
        do idx_prof=1,num_prof
          if (Time_profiles(idx_prof)>1) then
            sfl2l=sflow_parameters(idx_prof,6)
            sfl2r=sflow_parameters(idx_prof,7)
            sfl2c=sflow_parameters(idx_prof,8)
            k2sfl=sflow_parameters(idx_prof,9)
            v0sfl2=sflow_parameters(idx_prof,10)
            sfl3l=sflow_parameters(idx_prof,11)
            sfl3r=sflow_parameters(idx_prof,12)
            sfl3c=sflow_parameters(idx_prof,13)
            k3sfl=sflow_parameters(idx_prof,14)
            v0sfl3=sflow_parameters(idx_prof,15)
            tlsfl=sflow_parameters(idx_prof,16)
            trsfl=sflow_parameters(idx_prof,17)
            tcsfl=sflow_parameters(idx_prof,18)
            ktsfl=sflow_parameters(idx_prof,19)
            if (t+0.5*delta_t>=tlsfl .and. t+0.5*delta_t<=trsfl) then
              A_time=0.0
              if (Time_profiles(idx_prof)==2) then
                t_p=2.0*PI*ktsfl/(trsfl-tlsfl)
                A_time=0.5*(1.0-cos(t_p*(t+0.5*delta_t-tcsfl)))
              else if (Time_profiles(idx_prof)==3) then
                t_p=2.0*PI*ktsfl/(trsfl-tlsfl)
                A_time=sin(t_p*(t+0.5*delta_t-tcsfl))
              end if
              v_theta=v0sfl2*ARMS_Gauss(theta+0.5*delta_t*k_2_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))&
              &*ARMS_DGauss(phi+0.5*delta_t*k_2_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))&
              &/sin(PI*(theta+0.5*delta_t*k_2_th))/PI*A_time
              v_phi=v0sfl3*ARMS_DGauss(theta+0.5*delta_t*k_2_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))&
              &*ARMS_Gauss(phi+0.5*delta_t*k_2_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/PI*A_time
              k_3_th=k_3_th+v_theta/R/PI
              k_3_ph=k_3_ph+v_phi/R/PI/sin(PI*(theta+0.5*delta_t*k_2_th))
            end if
          end if
        end do
        do idx_prof=1,num_prof
          if (Time_profiles(idx_prof)>1) then
            sfl2l=sflow_parameters(idx_prof,6)
            sfl2r=sflow_parameters(idx_prof,7)
            sfl2c=sflow_parameters(idx_prof,8)
            k2sfl=sflow_parameters(idx_prof,9)
            v0sfl2=sflow_parameters(idx_prof,10)
            sfl3l=sflow_parameters(idx_prof,11)
            sfl3r=sflow_parameters(idx_prof,12)
            sfl3c=sflow_parameters(idx_prof,13)
            k3sfl=sflow_parameters(idx_prof,14)
            v0sfl3=sflow_parameters(idx_prof,15)
            tlsfl=sflow_parameters(idx_prof,16)
            trsfl=sflow_parameters(idx_prof,17)
            tcsfl=sflow_parameters(idx_prof,18)
            ktsfl=sflow_parameters(idx_prof,19)
            if (t+delta_t>=tlsfl .and. t+delta_t<=trsfl) then
              A_time=0.0
              if (Time_profiles(idx_prof)==2) then
                t_p=2.0*PI*ktsfl/(trsfl-tlsfl)
                A_time=0.5*(1.0-cos(t_p*(t+delta_t-tcsfl)))
              else if (Time_profiles(idx_prof)==3) then
                t_p=2.0*PI*ktsfl/(trsfl-tlsfl)
                A_time=sin(t_p*(t+delta_t-tcsfl))
              end if
              v_theta=v0sfl2*ARMS_Gauss(theta+delta_t*k_3_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))&
              &*ARMS_DGauss(phi+delta_t*k_3_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))&
              &/sin(PI*(theta+delta_t*k_3_th))/PI*A_time
              v_phi=v0sfl3*ARMS_DGauss(theta+delta_t*k_3_th,sfl2c,k2sfl/(sfl2r-sfl2l)/(sfl2r-sfl2l))&
              &*ARMS_Gauss(phi+delta_t*k_3_ph,sfl3c,k3sfl/(sfl3r-sfl3l)/(sfl3r-sfl3l))/PI*A_time
              k_4_th=k_4_th+v_theta/R/PI
              k_4_ph=k_4_ph+v_phi/R/PI/sin(PI*(theta+delta_t*k_3_th))
            end if
          end if
        end do
        
        arr_out(1)=theta+(k_1_th+2.0*k_2_th+2.0*k_3_th+k_4_th)*delta_t/6.0
        arr_out(2)=phi+(k_1_ph+2.0*k_2_ph+2.0*k_3_ph+k_4_ph)*delta_t/6.0
      end function surface_move_rk_4


      subroutine surfaceline_read_binary_cage(infilename,num_static_cage,len_static_cages,num_dynamic_cage,&
      &len_dynamic_cages,static_cages,dynamic_cages)
        character (len = 200) :: infilename
        integer(4) :: num_static_cage
        integer(4) :: num_dynamic_cage
        integer(4),allocatable :: len_static_cages(:)
        integer(4),allocatable :: len_dynamic_cages(:)
        real(4),allocatable :: static_cages(:,:)
        real(4),allocatable :: dynamic_cages(:,:)

        integer(4) :: idx,idx2,sum_static_cage,sum_dynamic_cage,static_cage_st,dynamic_cage_st
        sum_static_cage=0
        sum_dynamic_cage=0
        
        open (unit=1,file=infilename,form='unformatted',access='direct',recl=4)
        read (1,rec=1) num_static_cage
        allocate(len_static_cages(num_static_cage))
	do idx=1,num_static_cage
	  read (1,rec=idx+1) len_static_cages(idx)
	  sum_static_cage=sum_static_cage+len_static_cages(idx)
	end do
	
        read (1,rec=2+num_static_cage) num_dynamic_cage
        allocate(len_dynamic_cages(num_dynamic_cage))
	do idx=1,num_dynamic_cage
	  read (1,rec=idx+2+num_static_cage) len_dynamic_cages(idx)
	  sum_dynamic_cage=sum_dynamic_cage+len_dynamic_cages(idx)
	end do
	
	allocate(static_cages(2,sum_static_cage))
	allocate(dynamic_cages(2,sum_dynamic_cage))
	
	static_cage_st=0
	do idx=1,num_static_cage
	  do idx2=1,len_static_cages(idx)
	    read (1,rec=idx2+static_cage_st*2+2+num_static_cage+num_dynamic_cage) static_cages(1,idx2+static_cage_st)
	  end do
	  do idx2=1,len_static_cages(idx)
	    read (1,rec=idx2+static_cage_st*2+2+num_static_cage+num_dynamic_cage+len_static_cages(idx))&
	    & static_cages(2,idx2+static_cage_st)
	  end do
	  static_cage_st=static_cage_st+len_static_cages(idx)
	end do
	
	dynamic_cage_st=0
	do idx=1,num_dynamic_cage
	  do idx2=1,len_dynamic_cages(idx)
	    read (1,rec=idx2+dynamic_cage_st*2+2+num_static_cage+num_dynamic_cage+2*sum_static_cage)&
	    & dynamic_cages(1,idx2+dynamic_cage_st)
	  end do
	  do idx2=1,len_dynamic_cages(idx)
	    read (1,rec=idx2+dynamic_cage_st*2+2+num_static_cage+num_dynamic_cage+2*sum_static_cage+len_dynamic_cages(idx))&
	    & dynamic_cages(2,idx2+dynamic_cage_st)
	  end do
	  dynamic_cage_st=dynamic_cage_st+len_dynamic_cages(idx)
	end do
	
        close(1)
      end subroutine surfaceline_read_binary_cage


      subroutine surfaceline_write_binary_cage(outfilename,num_static_cage,len_static_cages,num_dynamic_cage,&
      &len_dynamic_cages,static_cages,dynamic_cages)
        character (len = 200) :: outfilename
        integer(4) :: num_static_cage
        integer(4) :: num_dynamic_cage
        integer(4) :: len_static_cages(:)
        integer(4) :: len_dynamic_cages(:)
        real(4) :: static_cages(:,:)
        real(4) :: dynamic_cages(:,:)

        integer(4) :: idx,idx2,sum_static_cage,static_cage_st,dynamic_cage_st
        sum_static_cage=0
        open (unit=2,file=outfilename,form='unformatted',status='REPLACE',recl=4)
	close(2)
        open (unit=2,file=outfilename,form='unformatted',access='direct',recl=4)
        write (2,rec=1) num_static_cage
	do idx=1,num_static_cage
	  write (2,rec=idx+1) len_static_cages(idx)
	  sum_static_cage=sum_static_cage+len_static_cages(idx)
	end do
	
        write (2,rec=2+num_static_cage) num_dynamic_cage
	do idx=1,num_dynamic_cage
	  write (2,rec=idx+2+num_static_cage) len_dynamic_cages(idx)
	end do
	
	static_cage_st=0
	do idx=1,num_static_cage
	  do idx2=1,len_static_cages(idx)
	    write (2,rec=idx2+static_cage_st*2+2+num_static_cage+num_dynamic_cage) static_cages(1,idx2+static_cage_st)
	  end do
	  do idx2=1,len_static_cages(idx)
	    write (2,rec=idx2+static_cage_st*2+2+num_static_cage+num_dynamic_cage+len_static_cages(idx))&
	    & static_cages(2,idx2+static_cage_st)
	  end do
	  static_cage_st=static_cage_st+len_static_cages(idx)
	end do
	
	dynamic_cage_st=0
	do idx=1,num_dynamic_cage
	  do idx2=1,len_dynamic_cages(idx)
	    write (2,rec=idx2+dynamic_cage_st*2+2+num_static_cage+num_dynamic_cage+2*sum_static_cage)&
	    & dynamic_cages(1,idx2+dynamic_cage_st)
	  end do
	  do idx2=1,len_dynamic_cages(idx)
	    write (2,rec=idx2+dynamic_cage_st*2+2+num_static_cage+num_dynamic_cage+2*sum_static_cage+len_dynamic_cages(idx))&
	    & dynamic_cages(2,idx2+dynamic_cage_st)
	  end do
	  dynamic_cage_st=dynamic_cage_st+len_dynamic_cages(idx)
	end do
	
        close(2)
      end subroutine surfaceline_write_binary_cage


      subroutine parse_SurfaceLine_CNT(CNT_file,ARMS_CNT_file,SurfaceLine_initial,t_start,delta_t,solar_Radius,&
      &num_t,plot_t,SurfaceLine_outputs)
        character (len = 200) :: CNT_file      !Control file for this program
        character (len = 200) :: ARMS_CNT_file !For the ARMS control file
        character (len = 200) :: SurfaceLine_initial
        real(4) :: t_start
        real(4) :: delta_t
        real(4) :: solar_Radius
        integer(4) :: num_t
        real(4),allocatable :: plot_t(:)
        character,allocatable :: SurfaceLine_outputs(:,:)

        character (len = 200) :: temp
        integer(4) :: idx,idx2

        open (unit=3,file=CNT_file,form="formatted")
        read (3,'(a)') temp
        read (3,'(a)') ARMS_CNT_file
        read (3,'(a)') temp
        read (3,'(a)') SurfaceLine_initial
        read (3,'(a)') temp
        read (3,*) t_start
        read (3,'(a)') temp
        read (3,*) delta_t
        read (3,'(a)') temp
        read (3,*) solar_Radius
        read (3,'(a)') temp
        read (3,*) num_t
	allocate(plot_t(num_t))
	allocate(SurfaceLine_outputs(200,num_t))
        read (3,'(a)') temp
	do idx=1,num_t
	  read (3,*) plot_t(idx)
	end do
        read (3,'(a)') temp
	do idx=1,num_t
	  read (3,'(a)') temp
	  do idx2=1,200
            SurfaceLine_outputs(idx2,idx)=temp(idx2:idx2)
	  end do
	end do
        close(3)
      end subroutine parse_SurfaceLine_CNT


      subroutine parse_sflow_CNT(ARMS_CNT_file,num_prof,Time_profiles,sflow_parameters)
        character (len = 200) :: ARMS_CNT_file
        integer(4) :: num_prof
        integer(4),allocatable :: Time_profiles(:)
        real(4),allocatable :: sflow_parameters(:,:)

        character (len = 200) :: line
        character (len = 200) :: split_line
        integer(4) :: num_lines,IO_stat,lines_read,idx_p
        logical :: reading_file, str_is_in
        
        num_lines=0
        reading_file=.true.
        open (unit=4,file=ARMS_CNT_file,form="formatted")
        do while (reading_file)
          read(4,*,iostat=IO_stat)
          if(IO_stat/=0) then 
            reading_file=.false.
          else
            num_lines=num_lines+1
          end if
        end do
        close(4)
        
        lines_read=0
        open (unit=4,file=ARMS_CNT_file,form="formatted")
        do while (lines_read<num_lines)
          read(4,'(a)') line
          lines_read=lines_read+1
          call str_lower(200,line)
          if (string_in_string(6,200,'nsflow',line)) then
            read(4,*) num_prof
            lines_read=lines_read+1
            allocate(Time_profiles(num_prof))
            allocate(sflow_parameters(num_prof,19))
            
            do idx_p=1,num_prof
              read(4,'(a)') line
              lines_read=lines_read+1
              read(4,*) sflow_parameters(idx_p,1), sflow_parameters(idx_p,2), &
              &sflow_parameters(idx_p,3), sflow_parameters(idx_p,4), sflow_parameters(idx_p,5)
              lines_read=lines_read+1
              read(4,'(a)') line
              lines_read=lines_read+1
              read(4,*) sflow_parameters(idx_p,6), sflow_parameters(idx_p,7), &
              &sflow_parameters(idx_p,8), sflow_parameters(idx_p,9), sflow_parameters(idx_p,10)
              lines_read=lines_read+1
              read(4,'(a)') line
              lines_read=lines_read+1
              read(4,*) sflow_parameters(idx_p,11), sflow_parameters(idx_p,12), &
              &sflow_parameters(idx_p,13), sflow_parameters(idx_p,14), sflow_parameters(idx_p,15)
              lines_read=lines_read+1
              read(4,'(a)') line
              lines_read=lines_read+1
              read(4,*) sflow_parameters(idx_p,16), sflow_parameters(idx_p,17), &
              &sflow_parameters(idx_p,18), sflow_parameters(idx_p,19)
              lines_read=lines_read+1
            end do
          end if
        end do
        close(4)
        
        idx_p=1
        lines_read=0
        open (unit=4,file=ARMS_CNT_file,form="formatted")
        do while (lines_read<num_lines)
          read(4,'(a)') line
          lines_read=lines_read+1
          call str_lower(200,line)
          if (string_in_string(23,200,'separable boundary flow',line)) then
            do while (lines_read<num_lines .and. idx_p<=num_prof)
              read(4,'(a)') line
              lines_read=lines_read+1
              if (string_stats_letter(line)) then
                call str_lower(200,line)
                call string_split_num(1,200,line,split_line)
                if (string_in_string(2,200,'ti',split_line)) then
		  call string_split_num(2,200,line,split_line)
                  if (split_line(1:1)=='u') then
                    Time_profiles(idx_p)=0
                  else if (split_line(1:1)=='l') then
                    Time_profiles(idx_p)=1
                  else if (split_line(1:1)=='c') then
                    Time_profiles(idx_p)=2
                  else if (split_line(1:1)=='s') then
                    Time_profiles(idx_p)=3
                  end if 
                end if
              else
                idx_p=idx_p+1
              end if
            end do
            
          end if
        end do
        close(4)
      end subroutine parse_sflow_CNT


      subroutine str_raise(len_str,string)
        character(*), intent(in out) :: string
        integer(4), intent(in) :: len_str
        integer(4) :: idx
      
        do idx=1,len_str
          select case(string(idx:idx))
            case("a":"z")
            string(idx:idx) = achar(iachar(string(idx:idx))-32)
          end select
        end do 
      end subroutine str_raise


      subroutine str_lower(len_str,string)
        character(*), intent(in out) :: string
        integer(4), intent(in) :: len_str
        integer(4) :: idx
 
        do idx=1,len_str
          select case(string(idx:idx))
            case("A":"Z")
              string(idx:idx) = achar(iachar(string(idx:idx))+32)
          end select
        end do  
      end subroutine str_lower


      !is_in = (string1 in string2)
      function string_in_string(len1,len2,string1,string2) result(is_in)
        integer(4), intent(in) :: len1
        integer(4), intent(in) :: len2
        character(*), intent(in) :: string1
        character(*), intent(in) :: string2
        logical :: is_in

        integer(4) :: idx,idx2
        character :: check_char1,check_char2
        logical :: is_in_temp
        
        if (len1>len2) then
          is_in=.false.
        else
          is_in=.false.
          idx=0
          do while (idx<=len2-len1)
            is_in_temp=.true.
            do idx2=1,len1
              check_char1=string1(idx2:idx2)
              check_char2=string2(idx2+idx:idx2+idx)
              if (check_char1/=check_char2) is_in_temp=.false.
            end do
            if (is_in_temp) then
              is_in=.true.
              idx=len2
            else
              idx=idx+1
            end if
          end do
        end if
      end function string_in_string


      function string_stats_letter(string) result(letter_start)
        character(*), intent(in) :: string
        logical :: letter_start

        if (iachar(string(1:1))>=65 .and. iachar(string(1:1))<=90) then
          letter_start=.true.
        else if (iachar(string(1:1))>=97 .and. iachar(string(1:1))<=122) then
          letter_start=.true.
        else
          letter_start=.false.
        end if
      end function string_stats_letter


      subroutine string_split_num(num,len_string,string_in,string_out)
        integer(4), intent(in) :: num
        integer(4), intent(in) :: len_string
        character(*), intent(in) :: string_in
        character(*) :: string_out
        
        integer(4) :: idx,idx2,num_current
        logical :: proceed
        
        do idx=1,len_string
          string_out(idx:idx)=' '
        end do
        
        idx=1
        proceed=.true.
        do while (idx<=len_string .and. proceed)
          if (iachar(string_in(idx:idx))>=65 .and. iachar(string_in(idx:idx))<=90) then
            proceed=.false.
          else if (iachar(string_in(idx:idx))>=97 .and. iachar(string_in(idx:idx))<=122) then
            proceed=.false.
          else
            idx=idx+1
          end if
        end do
        
        num_current=1
        do while (idx<=len_string .and. num_current<num)
          if (iachar(string_in(idx:idx))>=65 .and. iachar(string_in(idx:idx))<=90) then
            idx=idx+1
          else if (iachar(string_in(idx:idx))>=97 .and. iachar(string_in(idx:idx))<=122) then
            idx=idx+1
          else
            num_current=num_current+1
            proceed=.true.
            do while (idx<=len_string .and. proceed)
              if (iachar(string_in(idx:idx))>=65 .and. iachar(string_in(idx:idx))<=90) then
                proceed=.false.
              else if (iachar(string_in(idx:idx))>=97 .and. iachar(string_in(idx:idx))<=122) then
                proceed=.false.
              else
                idx=idx+1
              end if
            end do
          end if
        end do

        proceed=.true.
        idx2=1
        do while (idx<=len_string .and. proceed)
          if (iachar(string_in(idx:idx))>=65 .and. iachar(string_in(idx:idx))<=90) then
            string_out(idx2:idx2)=string_in(idx:idx)
            idx=idx+1
            idx2=idx2+1
          else if (iachar(string_in(idx:idx))>=97 .and. iachar(string_in(idx:idx))<=122) then
            string_out(idx2:idx2)=string_in(idx:idx)
            idx=idx+1
            idx2=idx2+1
          else
            proceed=.false.
          end if
        end do
      end subroutine string_split_num


      subroutine surfaceline_advance_cage(len_cage,cage_in,t_start,num_t,plot_t,delta_t,ARMS_CNT_file,&
      &solar_Radius,cage_out)
        integer(4), intent(in) :: len_cage
        real(4), intent(in) :: cage_in(:,:)
        real(4), intent(in) :: t_start
        integer(4), intent(in) :: num_t
        real(4) :: plot_t(:)
        real(4), intent(in) :: delta_t
        character (len = 200) :: ARMS_CNT_file
        real(4), intent(in) :: solar_Radius
        real(4),allocatable :: cage_out(:,:,:)
        
        real(4) :: t, t_end
        integer(4) :: idx_c,idx_t_out
        integer(4) :: num_prof
        integer(4),allocatable :: Time_profiles(:)
        real(4),allocatable :: sflow_parameters(:,:)
        real(4),allocatable :: cage(:,:)
        real(4) :: delta_t_current,new_pos(2)
        
        call parse_sflow_CNT(ARMS_CNT_file,num_prof,Time_profiles,sflow_parameters)
        t_end=plot_t(num_t)
        allocate(cage(2,len_cage))
        allocate(cage_out(2,len_cage,num_t))
        do idx_c=1,len_cage
          cage(1,idx_c)=1.0-cage_in(2,idx_c)/180.0
          cage(2,idx_c)=cage_in(1,idx_c)/180.0
        end do
        t=t_start
        idx_t_out=1
        
        do while (t<t_end)
          if (t<plot_t(idx_t_out)-delta_t) then
            delta_t_current=delta_t
            do idx_c=1,len_cage
              new_pos=surface_move_rk_4(solar_Radius,cage(1,idx_c),cage(2,idx_c),t,delta_t_current,&
              &num_prof,Time_profiles,sflow_parameters)
              cage(1,idx_c)=new_pos(1)
              cage(2,idx_c)=new_pos(2)
            end do
          else
            delta_t_current=plot_t(idx_t_out)-t
            do idx_c=1,len_cage
              new_pos=surface_move_rk_4(solar_Radius,cage(1,idx_c),cage(2,idx_c),t,delta_t_current,&
              &num_prof,Time_profiles,sflow_parameters)
              cage(1,idx_c)=new_pos(1)
              cage(2,idx_c)=new_pos(2)
              cage_out(1,idx_c,idx_t_out)=new_pos(2)*180.0
              cage_out(2,idx_c,idx_t_out)=(1.0-new_pos(1))*180.0
            end do
            idx_t_out=idx_t_out+1
          end if
          t=t+delta_t_current
        end do
        
        deallocate(cage)
      end subroutine surfaceline_advance_cage


      subroutine reconnected_advance_cage(len_cage,cage,t_start,t_end,delta_t,num_prof,Time_profiles,&
      &sflow_parameters,solar_Radius)
        integer(4), intent(in) :: len_cage
        real(4), intent(in out) :: cage(:,:)
        real(4), intent(in) :: t_start
        real(4), intent(in) :: t_end
        real(4), intent(in) :: delta_t
        integer(4) :: num_prof
        integer(4),allocatable :: Time_profiles(:)
        real(4),allocatable :: sflow_parameters(:,:)
        real(4), intent(in) :: solar_Radius
        
        real(4) :: t,delta_t_signed,delta_t_current,new_pos(2)
        integer(4) :: idx_c
        
        delta_t_signed=-abs(delta_t)
        t=t_end
        
        do while (t>t_start)
          if (t>t_start-delta_t_signed) then
            delta_t_current=delta_t_signed
          else
            delta_t_current=t_start-t
          end if
          do idx_c=1,len_cage
            new_pos=surface_move_rk_4(solar_Radius,cage(1,idx_c),cage(2,idx_c),t,delta_t_current,&
            &num_prof,Time_profiles,sflow_parameters)
            cage(1,idx_c)=new_pos(1)
            cage(2,idx_c)=new_pos(2)
          end do
          t=t+delta_t_current
        end do
        
      end subroutine reconnected_advance_cage


      subroutine parse_ConnectivityMap_CNT(CNT_file,ARMS_CNT_file,start_filename,end_filename,&
      &connection_filename,start_time,end_time,delta_t,solar_Radius)
        character (len = 200) :: CNT_file      !Control file for this program
        character (len = 200) :: ARMS_CNT_file !For the ARMS control file
        character (len = 200) :: start_filename
        character (len = 200) :: end_filename
        character (len = 200) :: connection_filename
        real(4) :: start_time
        real(4) :: end_time
        real(4) :: delta_t
        real(4) :: solar_Radius

        character (len = 200) :: temp

        open (unit=7,file=CNT_file,form="formatted")
        read (7,'(a)') temp
        read (7,'(a)') ARMS_CNT_file
        read (7,'(a)') temp
        read (7,'(a)') start_filename
        read (7,'(a)') temp
        read (7,'(a)') end_filename
        read (7,'(a)') temp
        read (7,'(a)') connection_filename
        read (7,'(a)') temp
        read (7,*) start_time
        read (7,'(a)') temp
        read (7,*) end_time
        read (7,'(a)') temp
        read (7,*) delta_t
        read (7,'(a)') temp
        read (7,*) solar_Radius
        close(7)
      end subroutine parse_ConnectivityMap_CNT


      subroutine parse_QSL_Rbinfile(filename,R_actual,num_t,num_p,theta,phi,Q)
        character (len = 200) :: filename
        real(4) :: R_actual
        integer(4) :: num_t
        integer(4) :: num_p
        real(4),allocatable :: theta(:,:)
        real(4),allocatable :: phi(:,:)
        real(4),allocatable :: Q(:,:)

        real(4),parameter :: PI=4.0*atan(1.0)
        integer(4) :: num_points,filesize,idx_t,idx_p,idx
        real(8) :: phi_1,phi_temp,R_temp,theta_temp
        real(4) :: Q_temp

        inquire(file=filename,size=filesize)
        num_points=filesize/28

        open (unit=8,file=filename,form='unformatted',access='direct',recl=28)
        read (8,rec=1) phi_1,theta_temp,R_temp,Q_temp
        num_t=1
        idx_t=2
        do while (idx_t<=num_points)
          read (8,rec=idx_t) phi_temp,theta_temp,R_temp,Q_temp
          if (phi_temp==phi_1) then
            num_t=num_t+1
            idx_t=idx_t+1
          else
            idx_t=num_points+1
          end if
        end do
        close(8)
	
	num_p=num_points/num_t
	allocate(theta(num_p,num_t))
	allocate(phi(num_p,num_t))
	allocate(Q(num_p,num_t))

        idx=1	
        open (unit=8,file=filename,form='unformatted',access='direct',recl=28)
        do idx_p=1,num_p
          do idx_t=1,num_t
            read (8,rec=idx) phi_temp,theta_temp,R_temp,Q_temp
            theta(idx_p,idx_t)=real(theta_temp,4)+0.5*PI
            phi(idx_p,idx_t)=real(phi_temp,4)-PI
            Q(idx_p,idx_t)=real(Q_temp,4)
            idx=idx+1
          end do
        end do
        close(8)
        R_actual=real(R_temp,4)

      end subroutine parse_QSL_Rbinfile


      function find_angle_idx(target_angle,num_arr,arr) result(idx)
        real(4) :: target_angle
        integer(4),intent(in) :: num_arr
        real(4),intent(in) :: arr(:)
        integer(4) :: idx
        
        integer(4) :: idx1
        idx=1
        do idx1=1,num_arr-1
        if (arr(idx1)<=target_angle .and. arr(idx1+1)>=target_angle) then
          if (target_angle>=0.5*(arr(idx1)+arr(idx1+1))) then
            idx=idx1+1
          else
            idx=idx1
          end if
          exit
        end if
        end do
      end function 


      subroutine categorize_reconnected(ARMS_CNT_file,start_filename,end_filename,start_time,&
      &end_time,delta_t,solar_Radius,num_t,num_p,connection_map)
        character (len = 200) :: ARMS_CNT_file
        character (len = 200) :: start_filename
        character (len = 200) :: end_filename
        real(4) :: start_time
        real(4) :: end_time
        real(4) :: delta_t
        real(4) :: solar_Radius
        integer(4) :: num_t
        integer(4) :: num_p
        integer(4),allocatable :: connection_map(:,:)

        real(4),parameter :: PI=4.0*atan(1.0)
        real(4),allocatable :: theta_st(:,:)
        real(4),allocatable :: phi_st(:,:)
        real(4),allocatable :: Q_st(:,:)
        real(4),allocatable :: theta_end(:,:)
        real(4),allocatable :: phi_end(:,:)
        real(4),allocatable :: Q_end(:,:)
        real(4) :: R_actual,sfl2l,sfl2r,sfl3l,sfl3r
        integer(4) :: idx_t,idx_p,idx_t_new,idx_p_new,idx,len_cage
        integer(4) :: num_prof
        integer(4),allocatable :: Time_profiles(:)
        real(4),allocatable :: sflow_parameters(:,:)
        integer(4),allocatable :: cage_locations(:,:)
        real(4),allocatable :: cage(:,:)
        real(4),allocatable :: cage_out(:,:,:)
        integer(4),allocatable :: cage_st_idx(:,:)
        logical :: inside_flow
        
        call parse_sflow_CNT(ARMS_CNT_file,num_prof,Time_profiles,sflow_parameters)
        
        call parse_QSL_Rbinfile(start_filename,R_actual,num_t,num_p,theta_st,phi_st,Q_st)
        call parse_QSL_Rbinfile(end_filename,R_actual,num_t,num_p,theta_end,phi_end,Q_end)
        deallocate(theta_end)
        deallocate(phi_end)
        allocate(connection_map(num_p,num_t))
        allocate(cage_locations(num_p,num_t))
        
        do idx_p=1,num_p
          do idx_t=1,num_t
            connection_map(idx_p,idx_t)=0
            if (Q_st(idx_p,idx_t)>0.0) then
              connection_map(idx_p,idx_t)=connection_map(idx_p,idx_t)+1
            end if
            if (Q_end(idx_p,idx_t)>0.0) then
              connection_map(idx_p,idx_t)=connection_map(idx_p,idx_t)+2
            end if
          end do
        end do
        
        len_cage=0
        do idx_p=1,num_p
          do idx_t=1,num_t
            inside_flow=.false.
            do idx=1,num_prof
              if (Time_profiles(idx)>0) then
                sfl2l=sflow_parameters(idx,6)
                sfl2r=sflow_parameters(idx,7)
                sfl3l=sflow_parameters(idx,11)
                sfl3r=sflow_parameters(idx,12)
                if (1.0-theta_st(1,idx_t)/PI>=sfl2l .and. 1.0-theta_st(1,idx_t)/PI<=sfl2r &
                &.and. phi_st(idx_p,1)/PI>=sfl3l .and. phi_st(idx_p,1)/PI<=sfl3r) then
                  inside_flow=.true.
                end if
              end if
            end do
            if (inside_flow) then
              cage_locations(idx_p,idx_t)=1
              len_cage=len_cage+1
            else
              cage_locations(idx_p,idx_t)=0
            end if
          end do
        end do
        
        allocate(cage(2,len_cage))
        allocate(cage_st_idx(2,len_cage))
        idx=1
        do idx_p=1,num_p
          do idx_t=1,num_t
            if (cage_locations(idx_p,idx_t)==1) then
              cage(1,idx)=1.0-theta_st(1,idx_t)/PI
              cage(2,idx)=phi_st(idx_p,1)/PI
              cage_st_idx(1,idx)=idx_t
              cage_st_idx(2,idx)=idx_p
              idx=idx+1
            end if
          end do
        end do
        
        call reconnected_advance_cage(len_cage,cage,start_time,end_time,delta_t,num_prof,Time_profiles,&
        &sflow_parameters,solar_Radius)
        
        do idx=1,len_cage
          idx_t=cage_st_idx(1,idx)
          idx_p=cage_st_idx(2,idx)
          idx_t_new=find_angle_idx(PI*(1.0-cage(1,idx)),num_t,theta_st(1,:))
          idx_p_new=find_angle_idx(PI*cage(2,idx),num_p,phi_st(:,1))
          if (Q_st(idx_p_new,idx_t_new)>0.0 .and. Q_end(idx_p,idx_t)<0.0) then
            connection_map(idx_p,idx_t)=connection_map(idx_p,idx_t)+4
          else if (Q_st(idx_p_new,idx_t_new)<0.0 .and. Q_end(idx_p,idx_t)>0.0) then
            connection_map(idx_p,idx_t)=connection_map(idx_p,idx_t)+4
          end if
        end do
        
        deallocate(cage_locations)
        deallocate(Q_end)
        deallocate(theta_st)
        deallocate(phi_st)
        deallocate(Q_st)
      end subroutine categorize_reconnected


      subroutine save_connection_map(ARMS_CNT_file,start_filename,end_filename,connection_filename,&
      &start_time,end_time,delta_t,solar_Radius)
        character (len = 200) :: ARMS_CNT_file
        character (len = 200) :: start_filename
        character (len = 200) :: end_filename
        character (len = 200) :: connection_filename
        real(4) :: start_time
        real(4) :: end_time
        real(4) :: delta_t
        real(4) :: solar_Radius

        integer(4) :: num_t
        integer(4) :: num_p
        integer(4),allocatable :: connection_map(:,:)
        integer(4) :: idx_t,idx_p,idx

        call categorize_reconnected(ARMS_CNT_file,start_filename,end_filename,start_time,&
        &end_time,delta_t,solar_Radius,num_t,num_p,connection_map)

        idx=3
        open (unit=9,file=connection_filename,form='unformatted',status='REPLACE',recl=4)
	close(9)
        open (unit=9,file=connection_filename,form='unformatted',access='direct',recl=4)
        write(9,rec=1) num_p
        write(9,rec=2) num_t
        do idx_p=1,num_p
          do idx_t=1,num_t
            write(9,rec=idx) connection_map(idx_p,idx_t)
            idx=idx+1
          end do
        end do
        close(9)

        deallocate(connection_map)
      end subroutine save_connection_map

end module ASoT_Functions_Fortran

