--By SMITHY, Joint Probabilistic Data Association Filter 2024
--Please see here for the implementation: https://steamcommunity.com/sharedfiles/filedetails/?id=3388742860&searchtext=

--To Do:
--Optimisation:
	--Iterative search
	--clustering optimisation, might be able to do seperate searches for clusters of tracks distant from each other

--Revise track initiation and deletion:
	--handle cases where tracks coalesce
	--Revise track initiation logic, make it more like track deletion

--Time-since-detection filter

--UI: 
--map zoom
--target info

--SHORTHAND NOTATIONS
m=math
s,tau=screen,m.pi*2
PB,PN,ON,IN,IB=property.getBool,property.getNumber,output.setNumber,input.getNumber,input.getBool

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
function mat_trans(A)
	local r = {}
	for i = 1, #A[1] do
		r[i] = {}
		for j = 1, #A do
			r[i][j] = A[j][i]
		end
	end
	return r
end
function mat_op(A, B, op)
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
function mat_inv(m) --3x3
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

function mat_eigen(matrix) --chatgpt wrote most of this
	local a, b = matrix[1][1], matrix[1][3]
	local c, d = matrix[3][1], matrix[3][3]

	local trace = a + d
	local sqrt_discriminant = math.sqrt(trace^2 - 4 * (a * d - b * c))
	local value = {(trace + sqrt_discriminant)/2,(trace - sqrt_discriminant)/2}
	local vector = {}
	for i = 1,2 do
		local v1, v2
		if b ~= 0 then
			v1, v2 = value[i] - d, b
		elseif c ~= 0 then
			v1, v2 = c, value[i] - a
		else
			v1, v2 = 1, 0
		end
		local magnitude = math.sqrt(v1^2 + v2^2)
		vector[i] = {v1 / magnitude, v2 / magnitude}
	end
	return {m.sqrt(value[1]*gate), m.sqrt(value[2]*gate)}, vector
end

function build(layer, used, prob) --Finds all feasible joint events. Starts at the first track (layer) and assigning it a possible radar return, then moves to the next layer and assigns another non-conflicting radar return. At the same time it finds the product of all the liklihoods to find the probability that a joint event may occur.
	if layer > #srt then
		for i = 1,nmeasure do
			if used[i]~=nil then
				posterior[used[i] ][i] = posterior[used[i] ][i] + prob  

			end
		end
		sum = sum + prob

	else

		build(layer + 1, used, prob*miss) --assign track to clutter

		for row = 1,#T[srt[layer] ].G do --#T[srt[layer] ].G ~= #returns, except that some returns may be outside a track's validation gate, so T[srt[layer] ].G is a list of returns that specific track srt[layer] can see. srt[] ~= #T, except some tracks may be deleted (T[2]=nil) so srt is just a table with no nil values. 
			local detection = T[srt[layer] ].G[row]
			if detection[2] ~= 0 and used[detection[1] ]==nil then
				used[detection[1] ] = srt[layer]
				build(layer + 1, used, prob*(detection[2] > 0 and detection[2] or miss))
				used[detection[1] ] = nil

			end
		end
	end
end



max_rng = PN("Max Range") --Used for UI only
min_rng = PN("Min Range")
max_tgt = PN("Max Targets")
gap = PN("Sweep Period (Ticks)")
cov = PN("Initial Variance")
gate = PN("Gate Size")

tws = PB("Output XYZ:")
trails = PB("Trails")
slew_reset = PB("View Offset") and max_rng/2 or 0
gates = PB("Draw Gates")
Rnoise = 3000 --Min measurement noise distance
miss = 2*10^-10--chance of clutter

--Initialize Kalman
h=gap/60

T={} --Stores all tracks

I={} --Identity matrix
Q={} --process noise
H={} --observation matrix, 3x3 since we have new measurements in the form of a 3D position
for i = 1,6 do
	I[i]={}
	Q[i]={}
	H[i]=i<4 and {} or nil --M 4
	a = (i-1)//3+1
	for j = 1,6 do
		b = (j-1)//3+1 --M
		I[i][j] = i==j and 1 or 0
		Q[i][j] = (i-1)%3==(j-1)%3 and h^(4-a-b)/2^(m.max(1-a,0)+m.max(1-b,0)) or 0
		if H[i]~= nil then
			H[i][j] = i==j and 1 or 0

		end
	end
end

h=1/60
Q = mat_op(Q,PN("Process Noise")) --990
F = mat_op(I,1)
for i = 4,6 do
	F[i-3][i]=h
	if i>6 then
		F[i-6][i]=0.5*h^2

	end
end
F_t = mat_trans(F)
H_t = mat_trans(H)

returns={} --table containing all new radar returns from the current sweep
list={} --used for trails
srt={} --same as T, except all nil values removed
launch={} --contains the track index for all tracks that have been selected

lsweep = 0

lpress=false
spress=false
lspress=false

--UI JAZZ
scl = max_rng/1000
lx=PN("X")*32 --change this so it's variable
lz=PN("Y")*32
hx=lx/2
hz=lz/2
slew = {0,slew_reset}

function onTick()
--UI MAP SLEW AND ZOOM
--	scl = m.max(IN(28),0.2)
	mscl=hz/(scl*1000)
	gx = lx/2-slew[1]*mscl
	gz = lz/2+slew[2]*mscl


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
	O_t = mat_trans(O)

--RADAR RETURN PROCESSING
	for i = 1,5 do
		local ind = 6+i*4
		local Z = {IN(ind-3),IN(ind-2)*tau,IN(ind-1)*tau}
		if Z[1]>min_rng and IN(ind) == 0 then
			table.insert(returns,mat_op(mat_mult(O_t,{{Z[1]*m.cos(Z[3])*m.sin(Z[2])},{m.sin(Z[3])*Z[1]},{Z[1]*m.cos(Z[3])*m.cos(Z[2])}}),mpos,0)) --converts from local polar to global euclidean

		end
	end

--JOINT PROBABILISTIC DATA ASSOCIATION FILTER
	for k,t in pairs(T) do
		--State Prediction
		t.X = mat_mult(F,t.X) --Predicts track position and velocity one tick in the future
		t.Pp = mat_mult(mat_mult(F,t.Pp),F_t)

	end

	sweep = (IN(29)+0.5)%1
	if sweep<lsweep or complete == false then --only activates once whenever the radar completes a full sweep, complete ==false isn't needed anymore (used to be used for iterative search functionality)
		srt = {} --same as T, except all nil values removed
		nmatch = {} --used to check if a return has been assigned to any tracks
		posterior = {} --Weight used in the state update weighted-average, posterior[2][3] would be the probability that track 2 is assigned to measurement 3
		nmeasure = #returns --not really needed anymore, but it used to be used in the iterative version since some returns may be added after the iterative search began
		for k,t in pairs(T) do
			table.insert(srt,k)
			--State Process Noise Prediction
			t.Pp = mat_op(t.Pp,Q,0)

			--Measurement Noise Matrix (ignoring azimuth,elevation variance, and just assuming distance variance in 3D)
			local t1 = {{t.X[1][1]-mpos[1][1]},{t.X[2][1]-mpos[2][1]},{t.X[3][1]-mpos[3][1]}}
			t1 = (m.max(m.sqrt(t1[1][1]^2+t1[2][1]^2+t1[3][1]^2),Rnoise)*0.02)^2/12 -----T
			R = {
			{t1,0,0},
			{0,t1,0},
			{0,0,t1}
			}

			--Kalman Gain
			t.pp = mat_op(mat_mult(mat_mult(H,t.Pp),H_t),R,0)
			pp_i,pp_d = mat_inv(t.pp)
			t.K = mat_mult(mat_mult(t.Pp,H_t),pp_i)

			--Likelihood of association based on mahalanobis distance
			t.V={}
			t.G={}
			posterior[k] = {} --initilise this for the event hyhpothesis search later
			for i = 1,#returns do
				posterior[k][i] = 0

				t.V[i] = mat_op(returns[i],mat_mult(H,t.X),1) --Innovation
				V1 = mat_mult(mat_mult(mat_trans(t.V[i]),pp_i),t.V[i])[1][1] --Number of standard deviations innovation is away
				if V1<gate then --Likelihood/validation gate, if innovation mahalanobis distance is less than x stdevs away
					table.insert(t.G,{i,m.exp(-0.5*V1)/m.sqrt(m.abs(pp_d)*tau^2)}) --add likelihood to table specific to that track
					nmatch[i] = false --if a return is never set to false we know it wasn't assigned to any tracks. 

				end
			end
		end


		--Kalman Filter State Update, and Event Hypothesis Search
		complete = nmeasure==0 or #srt==0
		if not complete then
			sum = 0
			build(1, {}, 1) --event hypothesis search

			complete = true
			for k,t in pairs(T) do
				t.beta = 0
				Z = mat_op(returns[1],0)
				for i = 1,#t.V do
					posterior[k][i] = posterior[k][i]*(sum~=0 and 1/sum or 0) --normalize posterior matrix
					t.beta = t.beta+posterior[k][i] -- 1 - t.beta would be the probability that the track is assigned to cluttere

					Z = mat_op(Z,mat_op(t.V[i],posterior[k][i]),0) --weighted innovation

				end

				--Covariance update
				t.Pp = mat_op(t.Pp,mat_op(mat_mult(mat_mult(t.K,t.pp),mat_trans(t.K)),t.beta),1)
				Zs = mat_mult(Z,mat_trans(Z))
				Ps = mat_op(Zs,0)
				for i = 1,#t.V do
					Ps = mat_op(Ps, mat_op(mat_op(mat_mult(t.V[i],mat_trans(t.V[i])),posterior[k][i]),Zs,1), 0)

				end
				t.Pp = mat_op(t.Pp, mat_op(mat_mult(mat_mult(t.K,Ps),mat_trans(t.K)),.09) ,1) --Correction term for uncertain association, change if .09 is set to 0 then covariance will just be a circle instead of an ellipse

				--State Update
				t.X = mat_op(t.X,mat_mult(t.K,Z),0)

				--UI trails
				if trails then
					table.insert(list,t.X)
					if #list>4*#srt then
						table.remove(list,1)
	
					end
				end
			end
		end

		if #srt~=0 then
			--Track Deletion
			delete = {}
			for i = 1,#srt do
				local t = srt[i]
				if T[t].beta~=nil then
					if T[t].beta<0.01 then --nil error here for some reason
						T[t].time = T[t].time ~= nil and T[t].time+1 or 0
						if T[t].time>3 then
							delete[i] = t 

						end
					else
						T[t].time = 0

					end
				end
			end
			for _,v in pairs(delete) do
				T[v] = nil

			end
			--Track Initiation 
			for i = 1,nmeasure do
				if nmatch[i]~=false and #srt<max_tgt then
					T[#T+1] = {X=returns[i],Pp=mat_op(I,cov)}
					count = count+1
					break

				end
			end
		else
			if nmeasure>0 then
				T[1] = {X=returns[1],Pp=mat_op(I,cov)}
				count=2

			end
		end
		returns = {}

	end

	--LAUNCH AND SELECTION LOGIC
	touchx = IN(30)
	touchy = IN(31)

	fspress = IN(32) == 1
	if spress then
		if m.abs(touchy-lz*.6-2)<6 and m.abs(m.abs(touchx-hx)-lx/16)<6 then
			scl = touchx<hx and scl+.1 or m.max(scl-0.1,0.1)
		elseif not lspress then
			L = mat_op(mat_mult(O_t,{{(touchx-gx)/mscl},{0},{(gz-touchy)/mscl}}),mpos,0)
			near = {0,0}
			for k,t in pairs(T) do
				local rng = m.sqrt((t.X[1][1]-L[1][1])^2+(t.X[3][1]-L[3][1])^2)
				if rng<scl*200 and (rng<near[2] or near[2]==0) then
					near = {k,rng,{t.X[1][1]-mpos[1][1],t.X[3][1]-mpos[3][1]}}
	
				end
			end
		end
	end

	press = IN(27) == 1
	if near~=nil then
		if near[1]~=0 and T[near[1]] ~=nil then
			L = mat_mult(O,{{T[near[1]].X[1][1]-mpos[1][1]},{0},{T[near[1]].X[3][1]-mpos[3][1]}})
			slew = {L[1][1],L[3][1]}
			if press and not lpress then
				count1 = count1==nil and 1 or count1+1
				if tws then
					table.insert(launch,{near[1],count1})
					T[near[1] ].launch = true

				else
					launch[1] = {near[1],count1}

				end
			end
		else 
			slew = {0,slew_reset}
			near[1] = 0

		end
	end
	
	lpress = press --I hate this
	lspress = spress
	spress = fspress

	--CYCLE TRACK OUTPUT
	if #launch>0 then
		out = out==nil and 1 or out
		out = out<=#launch and out or 1
		if T[launch[out][1] ]~= nil then
			send = T[launch[out][1] ]
			if #send.X>3 then
				for i = 1,6 do
					ON(i,send.X[i][1])
		
				end
				ON(7,launch[out][2])
			end
			out = out+1
		else
			table.remove(launch,out)
	
		end
	end

	lsweep = sweep

end
function drawElp(x, y, a, b, tilt) --merge this with drawCircle at some point
	local ang = 2*m.pi/8
	px, py = x + m.cos(tilt)*a, y+m.sin(tilt)*a
	for i = 0, 8 do
		px_, py_ = a * m.cos(ang * i), b * m.sin(ang * i)
		px_, py_ = x + m.cos(tilt)*px_+m.sin(tilt)*py_, y + m.sin(tilt)*px_-m.cos(tilt)*py_
		screen.drawLine(px_, py_, px, py)
		px, py = px_, py_

	end
end
function drawCircle(x, y, r, split)
	local step = 2*m.pi/24
	for i = 0,24 do
		if i%2 == 0 or not split then
			local i = step*i-hdg
			local step = step/2
			screen.drawLine(x+math.cos(i-step)*r,y+math.sin(i-step)*r,x+math.cos(i+step)*r,y+math.sin(i+step)*r)

		end
	end
end
comp = {"E","N","W","S"}
function onDraw()
	--Background
	hdg = m.atan(O[1][1],O[1][3])
	s.setColor(0,.074,3)
	s.drawRectF(0,0,lx,lz)

	s.setColor(.53,40,40)
	s.drawRectF(gx-1,gz-1,2,2)

	s.setColor(61,61,61)
	drawCircle(gx-.5,gz-.5,max_rng*mscl,false)	

	s.setColor(31,31,31)
	drawCircle(gx-.5,gz-.5,max_rng*mscl/2,true)

	for i=1,12 do
		local j = hdg+(i-1)*m.pi/6
		local size = max_rng*0.9*mscl
		if (i-1)%3~=0 then
			s.drawRectF(gx+m.sin(j)*size-0.5,gz+m.cos(j)*size-0.5,1,1)
		else
			s.drawText(gx-1.5+m.sin(j)*size,gz-1.5+m.cos(j)*size,comp[i//3+1])
		end
	end


	--Tracks
	for k,t in pairs(T) do
		s.setColor((100-10*k)%255,(21*k)%255,120/k)
		local L = mat_mult(O,mat_op({{t.X[1][1]},{t.X[2][1]},{t.X[3][1]}},mpos,1)) --anytime you see this I'm converting from global space to local x,y,z
		if t.launch~=nil then
			s.drawText(gx+L[1][1]*mscl-1,gz-L[3][1]*mscl-2,"X")
		else
			s.drawText(gx+L[1][1]*mscl-1,gz-L[3][1]*mscl-2,"*")
		end


		--Validation Gate
		if gates and t.pp ~= nil then
			FX = mat_mult(F,t.X)
			values, vectors = mat_eigen(t.pp)
			drawElp(gx+L[1][1]*mscl,gz-L[3][1]*mscl,values[1]*mscl,values[2]*mscl,m.atan(vectors[1][2],vectors[1][1])-hdg)

		end
	end
	if trails then
		--Trails
		s.setColor(233,144,20)
		for k,p in ipairs(list) do
			local L = mat_mult(O,mat_op({{p[1][1]},{p[2][1]},{p[3][1]}},mpos,1))
			s.drawText(gx+L[1][1]*mscl-1,gz-L[3][1]*mscl-3,".")
	
		end
		--Returns
		s.setColor(40,170,20)
		for k,p in ipairs(returns) do
			local L = mat_mult(O,mat_op({{p[1][1]},{p[2][1]},{p[3][1]}},mpos,1))
			s.drawText(gx+L[1][1]*mscl-1,gz-L[3][1]*mscl-3,"+")
	
		end
	end

	s.setColor(31,31,31)
	s.drawRectF(hx-lx/16-2,lz*.6,5,5) --for the love of god merge these
	s.drawRectF(hx+lx/16-2,lz*.6,5,5)

	s.setColor(61,61,61)
	s.drawText(hx-lx/16-1,lz*.6,"-")
	s.drawText(hx+lx/16-1,lz*.6,"+")

end
