!To compile, run the following line (without leading exclamation mark):
!gfortran ASoT_Functions_Fortran.F90 Connectivity_Map_Generate.F90 -o Connectivity_Map_Generate

!Filenames, directories, lines in control files are capped at 200 characters by len=200 magic number

program Connectivity_Map_Generate
      
      use ASoT_Functions_Fortran
      
      implicit none
      
      character (len = 200) :: CNT_file      !Control file for this program
      character (len = 200) :: ARMS_CNT_file !For the ARMS control file
      character (len = 200) :: start_filename
      character (len = 200) :: end_filename
      character (len = 200) :: connection_filename
      real(4) :: start_time
      real(4) :: end_time
      real(4) :: delta_t
      real(4) :: solar_Radius
      
      CNT_file='Connectivity_Map_Generate.cnt'
      call parse_ConnectivityMap_CNT(CNT_file,ARMS_CNT_file,start_filename,end_filename,&
      &connection_filename,start_time,end_time,delta_t,solar_Radius)
      
      call save_connection_map(ARMS_CNT_file,start_filename,end_filename,connection_filename,&
      &start_time,end_time,delta_t,solar_Radius)
      
end program Connectivity_Map_Generate


