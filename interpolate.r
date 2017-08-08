.mecdf.interpolate = function (Fh, Fst, nr, nc, x, u)
{	n = nrow (u)
	v = numeric (n)
	for (i in 1:n)
		v [i] = Fh (nr, nc, x, u [i,], Fst)
	v
}

#Fst ignored
.mecdf.step = function (nr, nc, x, u, Fst)
{	k = rep (TRUE, nr)
	for (j in 1:nc) k = k & (x [,j] <= u [j])
	sum (k) / nr
}

#todo: simplify...
#todo: remove interpolation over vertices with same values
.mecdf.continuous = function (nr, nc, x, u, Fst)
{	k = rep (NA, nc)
	a = b = p = q = numeric (nc)
	for (j in 1:nc)
	{	i = (x [,j] <= u [j])
		n1 = sum (i)
		if (n1 == 0) k [j] = -Inf
		else if (n1 == nr) k [j] = max (x [,j])
		else
		{	x1 = x [i, j]
			a [j] = max (x1)
			if (u [j] == a [j]) k [j] = a [j]
			else
			{	x2 = x [!i, j] 
				b [j] = min (x2)
			}
		}
	}
	if (any (!is.na (k) & k == -Inf) ) 0
	else if (all (!is.na (k) ) ) Fst (nr, nc, x, k)
	else
	{	p = (u - a) / (b - a)
		q = 1 - p
		nvert = 2^sum (is.na (k) )
		w = 1
		ust = NULL
		vst = numeric (nvert)
		for (j in nc:1)
		{	if (is.na (k [j]) )
			{	w = c (q [j] * w, p [j] * w)
				ust1 = cbind (a [j], ust)
				ust2 = cbind (b [j], ust)
				ust = rbind (ust1, ust2)
			}
			else ust = cbind (k [j], ust)
		}
		for (i in 1:nvert) vst [i] = Fst (nr, nc, x, ust [i,])
		sum (w * vst)
	}
}

.mecdf.vertex = function (nr, nc, x, u)
{	k = rep (TRUE, nr)
	for (j in 1:nc) k = k & (x [,j] <= u [j])
	(sum (k) - 1) / (nr - 1)
}

#assume x values are sorted
.uecdf.interpolate = function (Fh, nr, x, u)
{	n = length (u)
	v = numeric (n)
	for (i in 1:n)
		v [i] = Fh (nr, x, u [i])
	v
}

.uecdf.step = function (nr, x, u)
	sum (x <= u) / nr

.uecdf.continuous = function (nr, x, u)
{	k = sum (x <= u)
	a = (k - 1) / (nr - 1)
	b = k / (nr - 1)
	if (u == x [k]) a
	else
	{	p = (u - x [k]) / (x [k + 1] - x [k])
		q = 1 - p
		q * a + p * b
	}
}

