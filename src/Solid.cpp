# include "include/IncomNSSolver.h"
# include <cmath>

Solid::Solid(const array<double, 2> & Center, CD & Radius, const Mesh & mesh)
{
	center = Center; R = Radius;

	uFlag.initShape(mesh.get_Nxu(), mesh.get_Nyu(), mesh.get_Nzu());
	vFlag.initShape(mesh.get_Nxv(), mesh.get_Nyv(), mesh.get_Nzv());
	wFlag.initShape(mesh.get_Nxw(), mesh.get_Nyw(), mesh.get_Nzw());

	uFlag.setConstant(false); vFlag.setConstant(false); wFlag.setConstant(false); 

	determineFlags(mesh);
	determineIdx(mesh)
}


int Solid::determineFlags(const Mesh & mesh)
{
	auto &xu = mesh.get_xu(), &yu = mesh.get_yu(), &zu = mesh.get_zu(); 
	auto &xv = mesh.get_xv(), &yv = mesh.get_yv(), &zv = mesh.get_zv(); 
	auto &xw = mesh.get_xw(), &yw = mesh.get_yw(), &zw = mesh.get_zw(); 

	auto fu = [this, &xu, &yu] (CI &i, CI &j, CI &k) -> void
	{
		if ((xu[i]*xu[i] + yu[j]*yu[j]) < R*R) uFlag(i, j, k) = true;
	};

	auto fv = [this, &xv, &yv] (CI &i, CI &j, CI &k) -> void
	{
		if ((xv[i]*xv[i] + yv[j]*yv[j]) < R*R) vFlag(i, j, k) = true;
	};

	auto fw = [this, &xw, &yw] (CI &i, CI &j, CI &k) -> void
	{
		if ((xw[i]*xw[i] + yw[j]*yw[j]) < R*R) wFlag(i, j, k) = true;
	};

	tripleLoop(0, mesh.get_Nxu(), 0, mesh.get_Nyu(), 0, mesh.get_Nzu(), fu);
	tripleLoop(0, mesh.get_Nxv(), 0, mesh.get_Nyv(), 0, mesh.get_Nzv(), fv);
	tripleLoop(0, mesh.get_Nxw(), 0, mesh.get_Nyw(), 0, mesh.get_Nzw(), fw);

	return 0;
}


int Solid::determineIdx(const Mesh & mesh)
{
	const VD &xu = mesh.get_xu(), &yu = mesh.get_yu(), &zu = mesh.get_zu(); 
	const VD &xv = mesh.get_xv(), &yv = mesh.get_yv(), &zv = mesh.get_zv(); 
	const VD &xw = mesh.get_xw(), &yw = mesh.get_yw(), &zw = mesh.get_zw(); 

	const double &dx = mesh.get_dx(), &dy = mesh.get_dy(), &dz = mesh.get_dz(); 


	auto fu = [this, &xu, &yu, &dx, &dy] (CI &i, CI &j, CI &k) -> void
	{
		double dist;

		if ( ! uFlag(i, j, k) )
		{

			if ( uFlag(i, j+1, k) )
			{
				uIdxBC.push_back({i, j, k});
				uCorIdxBC.push_back({i, j-1, k});
				dist = sqrt(xu[i] * xu[i] + yu[j] * yu[j]) - R;
				uDist.push_back({dist, dist+dy});
			} 
			else if ( uFlag(i, j-1, k) )
			{
				uIdxBC.push_back({i, j, k});
				uCorIdxBC.push_back({i, j+1, k});
				dist = sqrt(xu[i] * xu[i] + yu[j] * yu[j]) - R;
				uDist.push_back({dist, dist+dy});
			} 
			else if ( uFlag(i+1, j, k) )
			{
				uIdxBC.push_back({i, j, k});
				uCorIdxBC.push_back({i-1, j, k});
				dist = sqrt(xu[i] * xu[i] + yu[j] * yu[j]) - R;
				uDist.push_back({dist, dist+dx});
			} 
			else if ( uFlag(i-1, j, k) )
			{
				uIdxBC.push_back({i, j, k});
				uCorIdxBC.push_back({i+1, j, k});
				dist = sqrt(xu[i] * xu[i] + yu[j] * yu[j]) - R;
				uDist.push_back({dist, dist+dx});
			}
		}
	};


	auto fv = [this, &xv, &yv, &dx, &dy] (CI &i, CI &j, CI &k) -> void
	{
		double dist;

		if ( ! vFlag(i, j, k) )
		{

			if ( vFlag(i, j+1, k) )
			{
				vIdxBC.push_back({i, j, k});
				vCorIdxBC.push_back({i, j-1, k});
				dist = sqrt(xv[i] * xv[i] + yv[j] * yv[j]) - R;
				vDist.push_back({dist, dist+dy});
			} 
			else if ( vFlag(i, j-1, k) )
			{
				vIdxBC.push_back({i, j, k});
				vCorIdxBC.push_back({i, j+1, k});
				dist = sqrt(xv[i] * xv[i] + yv[j] * yv[j]) - R;
				vDist.push_back({dist, dist+dy});
			} 
			else if ( vFlag(i+1, j, k) )
			{
				vIdxBC.push_back({i, j, k});
				vCorIdxBC.push_back({i-1, j, k});
				dist = sqrt(xv[i] * xv[i] + yv[j] * yv[j]) - R;
				vDist.push_back({dist, dist+dx});
			} 
			else if ( vFlag(i-1, j, k) )
			{
				vIdxBC.push_back({i, j, k});
				vCorIdxBC.push_back({i+1, j, k});
				dist = sqrt(xv[i] * xv[i] + yv[j] * yv[j]) - R;
				vDist.push_back({dist, dist+dx});
			}
		}
	};


	auto fw = [this, &xw, &yw, &dx, &dy] (CI &i, CI &j, CI &k) -> void
	{
		double dist;

		if ( ! wFlag(i, j, k) )
		{

			if ( wFlag(i, j+1, k) )
			{
				wIdxBC.push_back({i, j, k});
				wCorIdxBC.push_back({i, j-1, k});
				dist = sqrt(xw[i] * xw[i] + yw[j] * yw[j]) - R;
				wDist.push_back({dist, dist+dy});
			} 
			else if ( wFlag(i, j-1, k) )
			{
				wIdxBC.push_back({i, j, k});
				wCorIdxBC.push_back({i, j+1, k});
				dist = sqrt(xw[i] * xw[i] + yw[j] * yw[j]) - R;
				wDist.push_back({dist, dist+dy});
			} 
			else if ( wFlag(i+1, j, k) )
			{
				wIdxBC.push_back({i, j, k});
				wCorIdxBC.push_back({i-1, j, k});
				dist = sqrt(xw[i] * xw[i] + yw[j] * yw[j]) - R;
				wDist.push_back({dist, dist+dx});
			} 
			else if ( wFlag(i-1, j, k) )
			{
				wIdxBC.push_back({i, j, k});
				wCorIdxBC.push_back({i+1, j, k});
				dist = sqrt(xw[i] * xw[i] + yw[j] * yw[j]) - R;
				wDist.push_back({dist, dist+dx});
			}
		}
	};

	tripleLoop(0, mesh.get_Nxu(), 0, mesh.get_Nyu(), 0, mesh.get_Nzu(), fu);
	tripleLoop(0, mesh.get_Nxv(), 0, mesh.get_Nyv(), 0, mesh.get_Nzv(), fv);
	tripleLoop(0, mesh.get_Nxw(), 0, mesh.get_Nyw(), 0, mesh.get_Nzw(), fw);

			
	return 0;
}
