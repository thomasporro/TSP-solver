 cmake_minimum_required(VERSION 2.8)

 project(tsp)
 set(SRC
         src/main.c
         src/tsp.c
         src/utils.c
         src/read_input.c
         src/plot.c src/genetic_utils.c)

 add_executable(${PROJECT_NAME} ${SRC})

 include_directories(include)

 set(CPLEX /opt/ibm/ILOG/CPLEX_Studio201/cplex)
 set(CONCORDELIB /opt/concorde)

 target_link_libraries(${PROJECT_NAME} -lcplex -lm -lpthread -ldl)
 target_link_libraries(${PROJECT_NAME} -L${CPLEX}/lib/x86-64_linux/static_pic)
 target_link_libraries(${PROJECT_NAME} ${CONCORDELIB}/concorde.a)
 target_include_directories(${PROJECT_NAME} PRIVATE ${CPLEX}/include/ilcplex ${CONCORDELIB})
