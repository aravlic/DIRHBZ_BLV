FC = mpiifort  -qmkl -mcmodel=large 

OBJ =  dirhbz.o dirhbz_vapor.o 

run: $(OBJ) 
	$(FC) -o run $(OBJ)
