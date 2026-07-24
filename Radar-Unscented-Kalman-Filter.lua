--By SMITHY, Unscented Kalman Filter 2024
--SHORTHAND NOTATIONS
m=math
s,tau=screen,m.pi*2
IN,IB,ON,PN=input.getNumber,input.getBool,output.setNumber,property.getNumber

--MATRIX FUNCTIONS
function mat_mult(A,B)
	local r = {}
	for i = 1,#A do
		r[i] = {}
		for j = 1,#B[1] do 
			r[i][j] = 0
			for k = 1,#B do
				r[i][j] = r[i][j] + A[i][k]*B[k][j]
			end
		end
	end
	return r
end
function mat_trans(A) --transpose
	local r = {}
	for i = 1, #A[1] do
		r[i] = {}
		for j = 1, #A do
			r[i][j] = A[j][i]
		end
	end
	return r
end
function mat_op(A, B, op) --Addition,sub,scale
	local r = {}
	for i = 1, #A do
		r[i] = {}
		for j = 1, #A[1] do
			if op == 0 then
				r[i][j] = A[i][j] + B[i][j]
			elseif op == 1 then
				if #B[i]==1 then
					r[i][j] = A[i][j] - B[i][1]
				else
					r[i][j] = A[i][j] - B[i][j]
				end
			else
				r[i][j] = A[i][j] * B

			end
		end
	end
	return r
end
function mat_inv(m) --3x3 only
t1=m[2][2]*m[3][3]-m[3][2]*m[2][3]
t2=m[2][3]*m[3][1]-m[2][1]*m[3][3]
t3=m[2][1]*m[3][2]-m[3][1]*m[2][2]
d=m[1][1]*t1+m[1][2]*t2+m[1][3]*t3
id=1/d
return {
	{t1*id,(m[1][3]*m[3][2]-m[1][2]*m[3][3])*id,(m[1][2]*m[2][3]-m[1][3]*m[2][2])*id},
	{t2*id,(m[1][1]*m[3][3]-m[1][3]*m[3][1])*id,(m[2][1]*m[1][3]-m[1][1]*m[2][3])*id},
	{t3*id,(m[3][1]*m[1][2]-m[1][1]*m[3][2])*id,(m[1][1]*m[2][2]-m[2][1]*m[1][2])*id}
}, d
end

--INITILISE KALMAN
initial_var = PN("Intital Covariance") --this isn't that important tbh, the covariance will eventually converge no matter what this is initially set to (based on the measurement noise and process noise), but if set high the kalman filter will be more sensitive to start off
process = PN("Process Noise") --if you want to make the filter more responsive increase this, by doing so you increase the covariance matrix
dst = PN("Time Since Detection")+1 --this is used for calculating the process noise and the measurement variance, since tsd was already taken, it's called dst
h = (dst)/60 --time step
h2=h^2

dom = 2*(dst+1)*(dst+2) --the variance for a single uniformly distributed variable is given by range^2/12, but since we're taking the midpoint between the min and max value we have to recalculate variance based on (var(min) + var(max))/2
t1 = (tau*0.002)^2/dom
R = {
{1,0,0}, --range variance (needs to be updated in onTick)
{0,t1,0}, --az
{0,0,t1} --el
}



--Process Noise
I={} --identity matrix will be used later for initilising the covariance matrix
Q={}
H={} --kind of like an observation matrix thing for UKF
for i = 1,9 do --V: set to 6 for constant velocity (CV)
	I[i]={}
	Q[i]={}
	H[i]=i<4 and {} or nil
	a = (i-1)//3 --V: add + 1 for Q to be random variance in velocity 
	for j = 1,9 do
		b = (j-1)//3 --V: add + 1 for CV
		I[i][j] = i==j and 1 or 0
		Q[i][j] = (i-1)%3==(j-1)%3 and h^(4-a-b)/2^(m.max(1-a,0)+m.max(1-b,0)) or 0 --M

	end
end
Q = mat_op(Q,process) --process is normally the random variance in acceleration, unless the Q initialisation was changed above

--State Transition
F = mat_op(I,1)
for i = 4,9 do --V: set 9 to 6 for CV
	F[i-3][i]=h
	F[i-3][i+3]=h2*0.5 --V: remove this line for constant velocity case

end

--UKF Sigma Weights
N=9 --set to #X, i.e. for a constant velocity model this would be 6
num_points = 2*N+1 --number of sigma points to use in transform, +1 because we're also using the mean in the transform
kappa = 3-N --set to 3 for normal distributions, I think you should increase this for a uniform though
Nkappa=N+kappa

W = {{}}
Wd = {}
for i = 1, num_points do
	W[1][i] = 1/(2*Nkappa)

	Wd[i] = {}
	for j = 1,num_points do
		Wd[i][j] = i==j and W[1][i] or 0

	end
end
W[1][1] = kappa/Nkappa
Wd[1][1] = kappa/Nkappa

--INITIALISE UI
list={}
scl = 32
gain = 32*2
lx=9*32
lz=5*32
gx=lx/2
gz=lz/2


function onTick()
--PHYSICS SENSOR DATA
	mpos={{IN(1)},{IN(2)},{IN(3)}}

	rx,ry,rz=IN(4),IN(5),IN(6)
	cx,cy,cz=m.cos(rx),m.cos(ry),m.cos(rz)
	sx,sy,sz=m.sin(rx),m.sin(ry),m.sin(rz)

	O = { --Rotation Matrix
	{cy*cz,cy*sz,-sy},
	{-cx*sz+sx*sy*cz,cx*cz+sx*sy*sz,sx*cy},
	{sx*sz+cx*sy*cz,-sx*cz+cx*sy*sz,cx*cy}
	}

	tsd = IN(10)
	data = {IN(7),IN(8)*tau,IN(9)*tau}

	if tsd>0 and radar~=nil then
		--Time since detection midpoint average
		go = true
		for i = 1,3 do 
			radar[i] = {m.min(data[i],radar[i][1]),m.max(data[i],radar[i][2])}

		end
	else

		if go then --when valid return do kalman
			go = false

			--MEASUREMENT MATRIX / tsd average stuff
			Z={}
			for i = 1,3 do
				Z[i] = {(radar[i][1]+radar[i][2])/2}
			end

			world = mat_op(mat_mult(mat_trans(O_tsd),{{Z[1][1]*m.cos(Z[3][1])*m.sin(Z[2][1])},{m.sin(Z[3][1])*Z[1][1]},{Z[1][1]*m.cos(Z[3][1])*m.cos(Z[2][1])}}),mpos_tsd,0) --used for UI and state model initilisation

			--TRACK INITIATION
			if X==nil then
				Pp = mat_op(I,initial_var) --Initilise covariance matrix using identity 
				X = world --x,y,z guess position
				for i = 4,9 do
					X[i] = {0} --sets velocity and acceleration to 0
	
				end
			end
	
			--KALMAN FILTER
			--Select Sigma (sample) Points
			L = mat_op(Pp,0)
			for i = 1, #L do --Decompose to find sqrt(Nk*Pp). <-- Pp is just the variance so this is how you'd get the standard deviation (i.e. sqrt(2*Pp) will be 2 standard deviations away from the mean). I've decided to have the sigma points spread 3 stdev away from the mean here
				for j = 1, i do
					local sum = 0
					for k = 1, j - 1 do
						sum = sum + L[i][k] * L[j][k]
					end
		
					if i == j then
						L[i][j] = m.sqrt(m.abs(Pp[i][i]*Nkappa - sum))
					else
						L[i][j] = (Pp[i][j]*Nkappa - sum) / L[j][j]
					end
				end
			end
			
			V = mat_op(X,1) --this'll contain all the sigma points, the data in the first column is the mean X, second column is the vector for the first sigma point, ...
			for i =1,N do
				for j =1,N do
					V[j][i+1] = X[j][1] + L[j][i]
					V[j][N+i+1] = X[j][1] - L[j][i]

				end
			end
	
			--Propgate Points
			V = mat_mult(F,V) --update sigma points based on velocity and acceleration
			
			--Mean/Covariance
			X = mat_mult(V,mat_trans(W)) -- you probs don't even need to do this, but instead do mat_mult(F,X) the math works out the same. But you need to do this if your state transition is non-linear (i.e. constant turn model)
			D = mat_op(V,X,1) --error between state space estimate and each sigma point
			Pp = mat_op(mat_mult(mat_mult(D,Wd),mat_trans(D)),Q,0) --Predict covariance and add process noise
	
			--Transform to local
			for j=1,num_points do
				local L1 = mat_mult(O_tsd,mat_op({{V[1][j]},{V[2][j]},{V[3][j]}},mpos_tsd,1))
				local l = L1[1][1]^2+L1[3][1]^2
				H[1][j] = m.sqrt(l+L1[2][1]^2) --Range
				H[2][j] = m.atan(L1[1][1],L1[3][1]) --Azimuth
				H[3][j] = m.atan(L1[2][1],m.sqrt(l)) --Elevation
	
			end
			H_ = mat_mult(H,mat_trans(W)) --same thing above in the mean/covariance section, but this time necessary since it's undergone a non-linear transform
			H = mat_op(H,H_,1)


			--Measurement Noise
			R[1][1] = (H_[1][1]*0.02)^2/dom
	
			--Kalman Gain
			pp = mat_op(mat_mult(mat_mult(H,Wd),mat_trans(H)),R,0)
			pp_i,pp_d = mat_inv(pp)
	
			K = mat_mult(mat_mult(D,Wd),mat_trans(H))
			K = mat_mult(K,pp_i)
	
	
			--State Update
			X = mat_op(X,mat_mult(K,mat_op(Z,H_,1)),0)
			Pp = mat_op(Pp,mat_mult(mat_mult(K,pp),mat_trans(K)),1)--]]





	
			--UI STUFF, TRAILS
			table.insert(list,{{X[1][1],X[3][1]},{world[1][1],world[3][1]}})
			if #list>100 then
				table.remove(list,1)
		
			end
			--OUTPUT STATE SPACE
			for i = 1,9 do 
			output.setNumber(i,X[i][1])
			end
	
		else
			X=nil --if not tsd>0 last tick then set X to nil since no valid target
	
		end

		--TSD average stuff
		O_tsd={O[1],O[2],O[3]} --need to save rotation matrix when tsd starts to increase
		mpos_tsd = {mpos[1],mpos[2],mpos[3]} --I think I have to do deepcopies here
		radar = {}
		for i = 1,3 do
			radar[i] = {data[i],data[i]}

		end
	end


--UI MAP SLEW AND ZOOM
scl=IN(14)
mx=X~=nil and X[1][1] or mpos[1][1]
mz=X~=nil and X[3][1] or mpos[3][1]
ON(28,mpos[1][1])
ON(29,mpos[2][1])
ON(30,scl)
ON(31,mx)
ON(32,mz)

mscl=scl/(gz/lx)
aposx,aposz=map.screenToMap(mx,mz,mscl,lx,lz,lx-gx,lz-gz)
mx,mz=map.mapToScreen(aposx,aposz,mscl,lx,lz,mpos[1][1],mpos[3][1])
mscl=scl*1000/gz

end

function onDraw()
	--Current State Space
	if X~=nil then
		s.setColor((100-10)%255,(21)%255,120)
		local q=mat_op({{-mpos[1][1]+X[1][1]},{mpos[3][1]-X[3][1]}},1/mscl)
		s.drawText(mx+q[1][1],mz+q[2][1],"*")


	--Past State Spaces
		s.setColor(233,144,20)
		for i = 1,#list do
			local q=mat_op({{-mpos[1][1]+list[i][1][1]},{mpos[3][1]-list[i][1][2]}},1/mscl)
			s.drawText(mx+q[1][1],mz+q[2][1],".")
	
	
		end

	--Returns
		s.setColor(40,170,20)
		for i = 1,#list do
			local q=mat_op({{-mpos[1][1]+list[i][2][1]},{mpos[3][1]-list[i][2][2]}},1/mscl)
			s.drawText(mx+q[1][1],mz+q[2][1],"+")

		end
	end
end
