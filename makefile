FC = mpiifort  -qmkl -mcmodel=large 

OBJ =  dirhbz_bj.o dirhbz_vapor.o 

run: $(OBJ) 
	$(FC) -o run $(OBJ)
