
int NSSolverEuler::updateGhost()
{
	for(auto it=mesh.BCs.begin(); it!=mesh.BCs.end(); ++it)
	{
		switch (it->direction)
		{
			case -1:
				switch (it->uType)
				{
					case 0:

						for(int j=0; j<mesh.Ny; ++j){
							for(int k=0; k<mesh.Nz; ++k){
								u(it->uBCIdx, j, k) = u(uBCcorIdx);

							}
						}
						break;

					case 1:
						
						for(int j=0; j<mesh.Ny; ++j){
							for(int k=0; k<mesh.Nz; ++k){
								u(it->uBCIdx, j, k) = - u(uBCcorIdx) ;

							}
						}
		}
	}
	return 0;
}
