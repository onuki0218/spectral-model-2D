  T  ª   k820309    9          19.0        [úÄ^                                                                                                          
       ../src/mod_analysis.F90 MOD_ANALYSIS          @                                         
                            @                              
                                                           
                                                          
                                                                                                                       #         @                                                      #ASPECT    #ASPECT_SQRT 	   #ASPECT_SQRT_RE 
   #TIME                                                   
                                                 	     
                                                 
     
                 
                                      
      #         @                                                     #NK    #NL    #MU    #ASPECT                                                                 
       p         5 r    n                                           1p          & p         5 r    n                                      1  & p         5 r    n                                          1     5 r    n                                      1   5 r    n                                          1                                        
                                      
               @                                                   
                &                   &                                                                                                        
                &                   &                                           #         @                                      	               #NUM_THREADS    n                        0              Comp_set_num_threads                                                                       #         @                                      	               #DYNAMIC_THREADS    n                            3              Comp_set_dynamic                                                                   #         @                                      	               #NESTED    n                             6              Comp_set_nested                                                                  %         @                                                          n                         8              Comp_get_num_threads                    %         @                                                          n                         :              Comp_get_max_threads                    %         @                                                          n                         <              Comp_get_thread_num                    %         @                                                          n                         >              Comp_get_num_procs                    %         @                                                          n                         @              Comp_in_parallel                    %         @                                                          n                         B              Comp_in_final                    %         @                                                          n                         D              Comp_get_dynamic                    %         @                                                          n                         F              Comp_get_nested                    %         @                                                          n                         H              Comp_get_thread_limit                    #         @                                       	               #MAX_LEVELS !   n                         K              Comp_set_max_active_levels                                                           !           %         @                              "                            n                         M              Comp_get_max_active_levels                    %         @                              #                            n                         O              Comp_get_level                    %         @                              $                            n                         Q              Comp_get_active_level                    %         @                              %                           #LEVEL &   n                      S              Comp_get_ancestor_thread_num                                                              &           %         @                              '                           #LEVEL (   n                      U              Comp_get_team_size                                                              (           #         @                                 )     	               #KIND *   #CHUNK_SIZE +   n                       Y              Comp_set_schedule                                                             *                                                  +           #         @                                 ,     	               #KIND -   #CHUNK_SIZE .   n                       [              Comp_get_schedule                                                              -                                                    .            %         @                              /                            n                         ]              Comp_get_proc_bind                    %         @                              0                            n                         _              Comp_get_num_places                    %         @                              1                           #PLACE_NUM 2   n                          b              Comp_get_place_num_procs                                                          2           #         @                                 3     	               #PLACE_NUM 4   #IDS 5   n                         e              Comp_get_place_proc_ids                                                           4                  @                                5                         p          1     1                   %         @                              6                            n                         g              Comp_get_place_num                    %         @                              7                            n                         i              Comp_get_partition_num_places                    #         @                                 8     	               #PLACE_NUMS 9   n                         l              Comp_get_partition_place_nums                           @                                9                         p          1     1                   %         @                              :                     
       n                         n              Comp_get_wtime                    %         @                              ;                     
       n                         p              Comp_get_wtick                    %         @                              <                            n                         r              Comp_get_default_device                    #         @                                 =     	               #DEVICE_NUM >   n                         u              Comp_set_default_device                                                           >           %         @                              ?                            n                         w              Comp_get_num_devices                    %         @                              @                            n                         y              Comp_get_num_teams                    %         @                              A                            n                         {              Comp_get_team_num                    %         @                              B                            n                         }              Comp_get_cancellation                    %         @                              C                            n                                       Comp_is_initial_device                    %         @                              D                            n                                       Comp_get_initial_device                    %         @                              E                            n                                       Comp_get_device_num                    %         @                              F                           #KIND G   #DEVICE_NUM H   n                                     Comp_pause_resource                                                             G                                                  H           %         @                              I                           #KIND J   n                                     Comp_pause_resource_all                                                             J           %         @                              K                            n                                        Comp_get_supported_active_levels                    #         @                                 L     	               #EVENT M   n                                    Comp_fulfill_event                                                              M           #         @                                N     	               #SVAR O   n                                     Comp_init_lock                                                              O            #         @                                P     	               #SVAR Q   n                                     Comp_destroy_lock                                                              Q            #         @                                R     	               #SVAR S   n                                     Comp_set_lock                                                              S            #         @                                T     	               #SVAR U   n                                     Comp_unset_lock                                                              U            %         @                             V                           #SVAR W   n                                     Comp_test_lock                                                              W            #         @                                X     	               #NVAR Y   n                                     Comp_init_nest_lock                                                              Y            #         @                                Z     	               #NVAR [   n                                     Comp_destroy_nest_lock                                                              [            #         @                                \     	               #NVAR ]   n                                     Comp_set_nest_lock                                                              ]            #         @                                ^     	               #NVAR _   n                                     Comp_unset_nest_lock                                                              _            %         @                             `                           #NVAR a   n                       ¡              Comp_test_nest_lock                                                              a            %         @                              b                            n                         £              Comp_get_max_task_priority                    #         @                                 c     	               #SVAR d   #HINT e   n                             ¦              Comp_init_lock_with_hint                                                        d                                                   e           #         @                                 f     	               #NVAR g   #HINT h   n                             ¨              Comp_init_nest_lock_with_hint                                                        g                                                   h           %         @                              i                           #COMMAND j   #MODIFIER k   #ARG l   n                             ­              Comp_control_tool                                                       j                                                  k                                                   l            #         @                                 m     	               #ALLOCATOR n   n                          °              Comp_destroy_allocator                                                          n           #         @                                 o     	               #ALLOCATOR p   n                          ²              Comp_set_default_allocator                                                          p           %         @                              q                            n                         ´              Comp_get_default_allocator                    #         @                                 r     	               #SIZE s   n                       ·              Ckmp_set_stacksize                                                             s           #         @                                 t     	               #SIZE u   n                       ¹              Ckmp_set_stacksize_s                                                             u           #         @                                 v     	               #MSEC w   n                       ¼              Ckmp_set_blocktime                                                             w           #         @                                 x     	                n                         ¾              Ckmp_set_library_serial                    #         @                                 y     	                n                         À              Ckmp_set_library_turnaround                    #         @                                 z     	                n                         Â              Ckmp_set_library_throughput                    #         @                                 {     	               #LIBNUM |   n                             Å              Ckmp_set_library                                                       |           #         @                                 }     	               #STRING ~   n                             È              Ckmp_set_defaults                ,       @                                ~                         p          1     1                           %         @                                                          n                         Ê              Ckmp_get_stacksize                    %         @                                                          n                         Ì              Ckmp_get_stacksize_s                    %         @                                                          n                         Î              Ckmp_get_blocktime                    %         @                                                          n                         Ð              Ckmp_get_library                    #         @                                      	               #NUM    n                        Ó              Ckmp_set_disp_num_buffers                                                                       %         @                                                         #MASK    n                       Ö              Ckmp_set_affinity                                                                          %         @                                                         #MASK    n                       Ø              Ckmp_get_affinity                                                                          %         @                                                          n                         Ú              Ckmp_get_affinity_max_proc                    #         @                                      	               #MASK    n                       Ü              Ckmp_create_affinity_mask                                                                          #         @                                      	               #MASK    n                       Þ              Ckmp_destroy_affinity_mask                                                                          %         @                                                         #PROC    #MASK    n                             á              Ckmp_set_affinity_mask_proc                                                                                                                      %         @                                                         #PROC    #MASK    n                             ã              Ckmp_unset_affinity_mask_proc                                                                                                                      %         @                                                         #PROC    #MASK    n                             å              Ckmp_get_affinity_mask_proc                                                                                                                      %         @                                                         #SIZE    n                       ç              Ckmp_malloc                                                                        %         @                                                         #SIZE    #ALIGNMENT    n                        ê              Ckmp_aligned_malloc                                                                                                                         %         @                                                         #NELEM    #ELSIZE    n                          î              Ckmp_calloc                                                                                                                       %         @                                                         #PTR     #SIZE ¡   n                      ð              Ckmp_realloc                                                                                                                 ¡           #         @                                 ¢     	               #PTR £   n                 	       ò              Ckmp_free                                                            £           #         @                                 ¤     	                n                         ô              Ckmp_set_warnings_on                    #         @                                 ¥     	                n                         ö              Ckmp_set_warnings_off                    %         @                              ¦                           #CANCELKIND §   n                         ù              Ckmp_get_cancellation_status                                                           §           #         @                                   ¨                    #TIME ©             
                                 ©     
             -      fn#fn    Í   @   j   MOD_COMMON      @   J   MOD_FFTW '   M  @   J   MOD_GOVERNING_EQUATION      @   J   MISK    Í  @       NL+MOD_COMMON      @       NK+MOD_COMMON 8   M         EQN_ASPECT_RATIO+MOD_GOVERNING_EQUATION ?   Ð  @   a   EQN_ASPECT_RATIO%ASPECT+MOD_GOVERNING_EQUATION D     @   a   EQN_ASPECT_RATIO%ASPECT_SQRT+MOD_GOVERNING_EQUATION G   P  @   a   EQN_ASPECT_RATIO%ASPECT_SQRT_RE+MOD_GOVERNING_EQUATION =     @   a   EQN_ASPECT_RATIO%TIME+MOD_GOVERNING_EQUATION .   Ð  l       EQN_MU+MOD_GOVERNING_EQUATION 1   <    a   EQN_MU%MU+MOD_GOVERNING_EQUATION 5   M  @   a   EQN_MU%ASPECT+MOD_GOVERNING_EQUATION      ¤       E+MOD_COMMON    1  ¤       Q+MOD_COMMON ,   Õ  ±       OMP_SET_NUM_THREADS+OMP_LIB 8     @   a   OMP_SET_NUM_THREADS%NUM_THREADS+OMP_LIB (   Æ  ±       OMP_SET_DYNAMIC+OMP_LIB 8   w	  @   a   OMP_SET_DYNAMIC%DYNAMIC_THREADS+OMP_LIB '   ·	  §       OMP_SET_NESTED+OMP_LIB .   ^
  @   a   OMP_SET_NESTED%NESTED+OMP_LIB ,   
  ¨       OMP_GET_NUM_THREADS+OMP_LIB ,   F  ¨       OMP_GET_MAX_THREADS+OMP_LIB +   î  §       OMP_GET_THREAD_NUM+OMP_LIB *     ¦       OMP_GET_NUM_PROCS+OMP_LIB (   ;  ¤       OMP_IN_PARALLEL+OMP_LIB %   ß  ¡       OMP_IN_FINAL+OMP_LIB (     ¤       OMP_GET_DYNAMIC+OMP_LIB '   $  £       OMP_GET_NESTED+OMP_LIB -   Ç  ©       OMP_GET_THREAD_LIMIT+OMP_LIB 2   p  ¶       OMP_SET_MAX_ACTIVE_LEVELS+OMP_LIB =   &  @   a   OMP_SET_MAX_ACTIVE_LEVELS%MAX_LEVELS+OMP_LIB 2   f  ®       OMP_GET_MAX_ACTIVE_LEVELS+OMP_LIB &     ¢       OMP_GET_LEVEL+OMP_LIB -   ¶  ©       OMP_GET_ACTIVE_LEVEL+OMP_LIB 4   _  »       OMP_GET_ANCESTOR_THREAD_NUM+OMP_LIB :     @   a   OMP_GET_ANCESTOR_THREAD_NUM%LEVEL+OMP_LIB *   Z  ±       OMP_GET_TEAM_SIZE+OMP_LIB 0     @   a   OMP_GET_TEAM_SIZE%LEVEL+OMP_LIB )   K  ·       OMP_SET_SCHEDULE+OMP_LIB .     @   a   OMP_SET_SCHEDULE%KIND+OMP_LIB 4   B  @   a   OMP_SET_SCHEDULE%CHUNK_SIZE+OMP_LIB )     ·       OMP_GET_SCHEDULE+OMP_LIB .   9  @   a   OMP_GET_SCHEDULE%KIND+OMP_LIB 4   y  @   a   OMP_GET_SCHEDULE%CHUNK_SIZE+OMP_LIB *   ¹  ¦       OMP_GET_PROC_BIND+OMP_LIB +   _  §       OMP_GET_NUM_PLACES+OMP_LIB 0     »       OMP_GET_PLACE_NUM_PROCS+OMP_LIB :   Á  @   a   OMP_GET_PLACE_NUM_PROCS%PLACE_NUM+OMP_LIB /     »       OMP_GET_PLACE_PROC_IDS+OMP_LIB 9   ¼  @   a   OMP_GET_PLACE_PROC_IDS%PLACE_NUM+OMP_LIB 3   ü     a   OMP_GET_PLACE_PROC_IDS%IDS+OMP_LIB *     ¦       OMP_GET_PLACE_NUM+OMP_LIB 5   &  ±       OMP_GET_PARTITION_NUM_PLACES+OMP_LIB 5   ×  ¹       OMP_GET_PARTITION_PLACE_NUMS+OMP_LIB @        a   OMP_GET_PARTITION_PLACE_NUMS%PLACE_NUMS+OMP_LIB &     ¢       OMP_GET_WTIME+OMP_LIB &   ¶  ¢       OMP_GET_WTICK+OMP_LIB /   X  «       OMP_GET_DEFAULT_DEVICE+OMP_LIB /      ³       OMP_SET_DEFAULT_DEVICE+OMP_LIB :   ¶   @   a   OMP_SET_DEFAULT_DEVICE%DEVICE_NUM+OMP_LIB ,   ö   ¨       OMP_GET_NUM_DEVICES+OMP_LIB *   !  ¦       OMP_GET_NUM_TEAMS+OMP_LIB )   D"  ¥       OMP_GET_TEAM_NUM+OMP_LIB -   é"  ©       OMP_GET_CANCELLATION+OMP_LIB .   #  ª       OMP_IS_INITIAL_DEVICE+OMP_LIB /   <$  «       OMP_GET_INITIAL_DEVICE+OMP_LIB +   ç$  §       OMP_GET_DEVICE_NUM+OMP_LIB +   %  Á       OMP_PAUSE_RESOURCE+OMP_LIB 0   O&  @   a   OMP_PAUSE_RESOURCE%KIND+OMP_LIB 6   &  @   a   OMP_PAUSE_RESOURCE%DEVICE_NUM+OMP_LIB /   Ï&  µ       OMP_PAUSE_RESOURCE_ALL+OMP_LIB 4   '  @   a   OMP_PAUSE_RESOURCE_ALL%KIND+OMP_LIB 8   Ä'  ´       OMP_GET_SUPPORTED_ACTIVE_LEVELS+OMP_LIB *   x(  ©       OMP_FULFILL_EVENT+OMP_LIB 0   !)  @   a   OMP_FULFILL_EVENT%EVENT+OMP_LIB &   a)  ¤       OMP_INIT_LOCK+OMP_LIB +   *  @   a   OMP_INIT_LOCK%SVAR+OMP_LIB )   E*  §       OMP_DESTROY_LOCK+OMP_LIB .   ì*  @   a   OMP_DESTROY_LOCK%SVAR+OMP_LIB %   ,+  £       OMP_SET_LOCK+OMP_LIB *   Ï+  @   a   OMP_SET_LOCK%SVAR+OMP_LIB '   ,  ¥       OMP_UNSET_LOCK+OMP_LIB ,   ´,  @   a   OMP_UNSET_LOCK%SVAR+OMP_LIB &   ô,  ¬       OMP_TEST_LOCK+OMP_LIB +    -  @   a   OMP_TEST_LOCK%SVAR+OMP_LIB +   à-  ©       OMP_INIT_NEST_LOCK+OMP_LIB 0   .  @   a   OMP_INIT_NEST_LOCK%NVAR+OMP_LIB .   É.  ¬       OMP_DESTROY_NEST_LOCK+OMP_LIB 3   u/  @   a   OMP_DESTROY_NEST_LOCK%NVAR+OMP_LIB *   µ/  ¨       OMP_SET_NEST_LOCK+OMP_LIB /   ]0  @   a   OMP_SET_NEST_LOCK%NVAR+OMP_LIB ,   0  ª       OMP_UNSET_NEST_LOCK+OMP_LIB 1   G1  @   a   OMP_UNSET_NEST_LOCK%NVAR+OMP_LIB +   1  ±       OMP_TEST_NEST_LOCK+OMP_LIB 0   82  @   a   OMP_TEST_NEST_LOCK%NVAR+OMP_LIB 2   x2  ®       OMP_GET_MAX_TASK_PRIORITY+OMP_LIB 0   &3  ¸       OMP_INIT_LOCK_WITH_HINT+OMP_LIB 5   Þ3  @   a   OMP_INIT_LOCK_WITH_HINT%SVAR+OMP_LIB 5   4  @   a   OMP_INIT_LOCK_WITH_HINT%HINT+OMP_LIB 5   ^4  ½       OMP_INIT_NEST_LOCK_WITH_HINT+OMP_LIB :   5  @   a   OMP_INIT_NEST_LOCK_WITH_HINT%NVAR+OMP_LIB :   [5  @   a   OMP_INIT_NEST_LOCK_WITH_HINT%HINT+OMP_LIB )   5  É       OMP_CONTROL_TOOL+OMP_LIB 1   d6  @   a   OMP_CONTROL_TOOL%COMMAND+OMP_LIB 2   ¤6  @   a   OMP_CONTROL_TOOL%MODIFIER+OMP_LIB -   ä6  @   a   OMP_CONTROL_TOOL%ARG+OMP_LIB .   $7  ±       OMP_DESTROY_ALLOCATOR+OMP_LIB 8   Õ7  @   a   OMP_DESTROY_ALLOCATOR%ALLOCATOR+OMP_LIB 2   8  µ       OMP_SET_DEFAULT_ALLOCATOR+OMP_LIB <   Ê8  @   a   OMP_SET_DEFAULT_ALLOCATOR%ALLOCATOR+OMP_LIB 2   
9  ®       OMP_GET_DEFAULT_ALLOCATOR+OMP_LIB *   ¸9  ¨       KMP_SET_STACKSIZE+OMP_LIB /   `:  @   a   KMP_SET_STACKSIZE%SIZE+OMP_LIB ,    :  ª       KMP_SET_STACKSIZE_S+OMP_LIB 1   J;  @   a   KMP_SET_STACKSIZE_S%SIZE+OMP_LIB *   ;  ¨       KMP_SET_BLOCKTIME+OMP_LIB /   2<  @   a   KMP_SET_BLOCKTIME%MSEC+OMP_LIB /   r<  £       KMP_SET_LIBRARY_SERIAL+OMP_LIB 3   =  §       KMP_SET_LIBRARY_TURNAROUND+OMP_LIB 3   ¼=  §       KMP_SET_LIBRARY_THROUGHPUT+OMP_LIB (   c>  ¨       KMP_SET_LIBRARY+OMP_LIB /   ?  @   a   KMP_SET_LIBRARY%LIBNUM+OMP_LIB )   K?  ©       KMP_SET_DEFAULTS+OMP_LIB 0   ô?     a   KMP_SET_DEFAULTS%STRING+OMP_LIB *   @  ¦       KMP_GET_STACKSIZE+OMP_LIB ,   &A  ¨       KMP_GET_STACKSIZE_S+OMP_LIB *   ÎA  ¦       KMP_GET_BLOCKTIME+OMP_LIB (   tB  ¤       KMP_GET_LIBRARY+OMP_LIB 1   C  ®       KMP_SET_DISP_NUM_BUFFERS+OMP_LIB 5   ÆC  @   a   KMP_SET_DISP_NUM_BUFFERS%NUM+OMP_LIB )   D  ¯       KMP_SET_AFFINITY+OMP_LIB .   µD  @   a   KMP_SET_AFFINITY%MASK+OMP_LIB )   õD  ¯       KMP_GET_AFFINITY+OMP_LIB .   ¤E  @   a   KMP_GET_AFFINITY%MASK+OMP_LIB 2   äE  ®       KMP_GET_AFFINITY_MAX_PROC+OMP_LIB 1   F  ¯       KMP_CREATE_AFFINITY_MASK+OMP_LIB 6   AG  @   a   KMP_CREATE_AFFINITY_MASK%MASK+OMP_LIB 2   G  °       KMP_DESTROY_AFFINITY_MASK+OMP_LIB 7   1H  @   a   KMP_DESTROY_AFFINITY_MASK%MASK+OMP_LIB 3   qH  Ã       KMP_SET_AFFINITY_MASK_PROC+OMP_LIB 8   4I  @   a   KMP_SET_AFFINITY_MASK_PROC%PROC+OMP_LIB 8   tI  @   a   KMP_SET_AFFINITY_MASK_PROC%MASK+OMP_LIB 5   ´I  Å       KMP_UNSET_AFFINITY_MASK_PROC+OMP_LIB :   yJ  @   a   KMP_UNSET_AFFINITY_MASK_PROC%PROC+OMP_LIB :   ¹J  @   a   KMP_UNSET_AFFINITY_MASK_PROC%MASK+OMP_LIB 3   ùJ  Ã       KMP_GET_AFFINITY_MASK_PROC+OMP_LIB 8   ¼K  @   a   KMP_GET_AFFINITY_MASK_PROC%PROC+OMP_LIB 8   üK  @   a   KMP_GET_AFFINITY_MASK_PROC%MASK+OMP_LIB #   <L  ©       KMP_MALLOC+OMP_LIB (   åL  @   a   KMP_MALLOC%SIZE+OMP_LIB +   %M  À       KMP_ALIGNED_MALLOC+OMP_LIB 0   åM  @   a   KMP_ALIGNED_MALLOC%SIZE+OMP_LIB 5   %N  @   a   KMP_ALIGNED_MALLOC%ALIGNMENT+OMP_LIB #   eN  ¶       KMP_CALLOC+OMP_LIB )   O  @   a   KMP_CALLOC%NELEM+OMP_LIB *   [O  @   a   KMP_CALLOC%ELSIZE+OMP_LIB $   O  ³       KMP_REALLOC+OMP_LIB (   NP  @   a   KMP_REALLOC%PTR+OMP_LIB )   P  @   a   KMP_REALLOC%SIZE+OMP_LIB !   ÎP         KMP_FREE+OMP_LIB %   lQ  @   a   KMP_FREE%PTR+OMP_LIB ,   ¬Q          KMP_SET_WARNINGS_ON+OMP_LIB -   LR  ¡       KMP_SET_WARNINGS_OFF+OMP_LIB 4   íR  À       KMP_GET_CANCELLATION_STATUS+OMP_LIB ?   ­S  @   a   KMP_GET_CANCELLATION_STATUS%CANCELKIND+OMP_LIB $   íS  R       ANALYSIS_CAL_ENERGY )   ?T  @   a   ANALYSIS_CAL_ENERGY%TIME 