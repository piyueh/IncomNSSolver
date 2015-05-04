# include "include/IncomNSSolver.h"
# include <cmath>

Solid::Solid(const array<double, 2> & Center, CD & Radius, const Mesh & mesh)
{
	center = Center; R = Radius;

	uFlag.initShape(-1, mesh.get_Nxu(), -1, mesh.get_Nyu(), -1, mesh.get_Nzu());
	vFlag.initShape(-1, mesh.get_Nxv(), -1, mesh.get_Nyv(), -1, mesh.get_Nzv());
	wFlag.initShape(-1, mesh.get_Nxw(), -1, mesh.get_Nyw(), -1, mesh.get_Nzw());

	uFlag.setConstant(false); vFlag.setConstant(false); wFlag.setConstant(false); 

	determineFlags(mesh);
	determineIdx(mesh);
}


int Solid::determineFlags(const Mesh & mesh)
{
	VD xu = mesh.get_xu(), yu = mesh.get_yu(), zu = mesh.get_zu(); 
	VD xv = mesh.get_xv(), yv = mesh.get_yv(), zv = mesh.get_zv(); 
	VD xw = mesh.get_xw(), yw = mesh.get_yw(), zw = mesh.get_zw(); 

	auto fx = [this] (double & d) -> void { d -= center[0]; };
	auto fy = [this] (double & d) -> void { d -= center[1]; };

	for_each(xu.begin(), xu.end(), fx);
	for_each(xv.begin(), xv.end(), fx);
	for_each(xw.begin(), xw.end(), fx);

	for_each(yu.begin(), yu.end(), fy);
	for_each(yv.begin(), yv.end(), fy);
	for_each(yw.begin(), yw.end(), fy);

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
	VD xu = mesh.get_xu(), yu = mesh.get_yu(), zu = mesh.get_zu(); 
	VD xv = mesh.get_xv(), yv = mesh.get_yv(), zv = mesh.get_zv(); 
	VD xw = mesh.get_xw(), yw = mesh.get_yw(), zw = mesh.get_zw(); 

	auto fx = [this] (double & d) -> void { d -= center[0]; };
	auto fy = [this] (double & d) -> void { d -= center[1]; };


	for_each(xu.begin(), xu.end(), fx);
	for_each(xv.begin(), xv.end(), fx);
	for_each(xw.begin(), xw.end(), fx);

	for_each(yu.begin(), yu.end(), fy);
	for_each(yv.begin(), yv.end(), fy);
	for_each(yw.begin(), yw.end(), fy);

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

	NPu = uIdxBC.size();
	NPv = vIdxBC.size();
	NPw = wIdxBC.size();

			
	return 0;
}


int Solid::updVelocity(A3Dd &u, A3Dd &v, A3Dd &w)
{
	int idx1, idx2, idx3;
	int cIdx1, cIdx2, cIdx3;

	double uBase;

	for(int i=0; i<NPu; ++i)
	{
		idx1 = uIdxBC[i][0];
		idx2 = uIdxBC[i][1];
		idx3 = uIdxBC[i][2];

		cIdx1 = uCorIdxBC[i][0];
		cIdx2 = uCorIdxBC[i][1];
		cIdx3 = uCorIdxBC[i][2];

		uBase = u(cIdx1, cIdx2, cIdx3);

		u(idx1, idx2, idx3) = uBase * uDist[i][0] / uDist[i][1];
	}

	for(int i=0; i<NPv; ++i)
	{
		idx1 = vIdxBC[i][0];
		idx2 = vIdxBC[i][1];
		idx3 = vIdxBC[i][2];

		cIdx1 = vCorIdxBC[i][0];
		cIdx2 = vCorIdxBC[i][1];
		cIdx3 = vCorIdxBC[i][2];

		uBase = v(cIdx1, cIdx2, cIdx3);

		v(idx1, idx2, idx3) = uBase * vDist[i][0] / vDist[i][1];
	}

	for(int i=0; i<NPw; ++i)
	{
		idx1 = wIdxBC[i][0];
		idx2 = wIdxBC[i][1];
		idx3 = wIdxBC[i][2];

		cIdx1 = wCorIdxBC[i][0];
		cIdx2 = wCorIdxBC[i][1];
		cIdx3 = wCorIdxBC[i][2];

		uBase = w(cIdx1, cIdx2, cIdx3);

		w(idx1, idx2, idx3) = uBase * wDist[i][0] / wDist[i][1];
	}

	return 0;
}
