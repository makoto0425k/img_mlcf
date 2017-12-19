source("parallel.R")

args <- commandArgs(trailingOnly = T)

filename <- args[1]
cf.no <- as.numeric(args[2])
inner.no <- NA
layer <- 4

if (length(args) == 3) {
	inner.no <- as.numeric(args[3])
	layer <- 5
}

debug <- FALSE

if (length(args) == 0) {
	source("defaultargs.R")
}

# ポータル一覧の読み込み

list.filename <- paste(filename, ".txt", sep="")
if (!file.exists(list.filename)) {
	list.filename <- paste(filename, ".txt.done", sep="")
	if (!file.exists(list.filename)) {
    quit()
  }
}

portal.list <- read.csv(list.filename, header=F, col.names=c("name","lat","lng"), stringsAsFactors=F)
portals <- cbind(cf.no = 1:nrow(portal.list), as.matrix(portal.list[, c("lat","lng")]))
portal.names <- portal.list[, "name"]
rm(portal.list)

# 見つかった均一多重の構成の読み込み

if (layer == 4) {
	ext <- ".mlcf4"
	proc.filename <- paste(filename, ".", cf.no, ext, ".path.txt", sep="")
}
if (layer == 5) {
	ext <- ".mlcf5"
	proc.filename <- paste(filename, ".", cf.no, ".", inner.no, ext, ".path.txt", sep="")
}

ext <- paste("mlcf", layer, sep="")
if (!file.exists(paste(filename, ".", ext, ".rdat", sep=""))) {
	quit()
}

load(paste(filename, ".", ext, ".rdat", sep=""))


# 対象の均一多重の選択

mlcf <- mlcfs[[cf.no]]

if (layer == 5) {
  # 均一5重の場合は、内部の均一4重の一つを選択
	vertexes <- matrix(c(1,2, 2,3, 3,1), byrow=T, ncol=2)

	inner.vertexes <- c(mlcfs[[cf.no]]$vertexes[vertexes[inner.no, ]], mlcf$inner[[1]]$center)
	inner.inner <- c(mlcf$inner[[1]]$list1, mlcf$inner[[1]]$list2, mlcf$inner[[1]]$list3)[inner.no]
	mlcf <- list(vertexes = inner.vertexes, inner = inner.inner)
}

inner <- mlcf$inner[[1]]


# 構成番号とポータル番号を対応させておく

mlcf.portals <- c(
	mlcf$vertexes[1],
	mlcf$vertexes[2],
	mlcf$vertexes[3],
	inner$center,
	inner$list1[[1]]$center,
	inner$list1[[1]]$list1[1],
	inner$list1[[1]]$list2[1],
	inner$list1[[1]]$list3[1],
	inner$list2[[1]]$center,
	inner$list2[[1]]$list1[1],
	inner$list2[[1]]$list2[1],
	inner$list2[[1]]$list3[1],
	inner$list3[[1]]$center,
	inner$list3[[1]]$list1[1],
	inner$list3[[1]]$list2[1],
	inner$list3[[1]]$list3[1]
) 

mlcf.posid <- c(
	1,2,3,9,91,911,912,913,92,921,922,923,93,931,932,933
)


# 均一多重を構成するポータルの座標情報

lnglat <- c()
for (p1 in mlcf.portals) {
	lnglat <- c(lnglat, portals[p1, c("lng", "lat")])
}
lnglat <- matrix(lnglat, byrow=T, ncol=2)
colnames(lnglat) <- c("lng", "lat")


# 距離マトリクス

source("lnglatdist.R")

portal.dists <- c()
for (i in 1:nrow(lnglat)) {
	portal.dists <- c(portal.dists, lnglat.dist(lnglat[i, "lng"], lnglat[i, "lat"], lnglat[, "lng"], lnglat[, "lat"]))
}
portal.dists <- matrix(portal.dists, byrow=T, nrow = nrow(lnglat))
rownames(portal.dists) <- portal.names[mlcf.portals]



if (file.exists(proc.filename)) {

	proc.file <- scan(proc.filename, what = character(), sep = "\n", blank.lines.skip = F)
	
	blank <- proc.file == ""
	blank2 <- which(c(blank) & c(FALSE, blank[-length(blank)]) )
	
	path.patterns <- list()
	path <- c()
	for (pattern.index in 1:6) { 
		line <- blank2[1] + 1 + 17 * (pattern.index - 1)
		proc.path <- unlist(strsplit(proc.file[line:(line+15)], "\t"))
		proc.path <- proc.path[proc.path != "<NG>"]
		dim(proc.path) <- c(2,16)
		path.patterns[[pattern.index]] <- data.frame(posid = as.numeric(proc.path[1, ]), name = proc.path[2, ], stringsAsFactors = F)
	}

} else {
	
	# クラスタ分析
	
	#plot(hclust(as.dist(portal.dists), method="centroid"))
	#plot(hclust(as.dist(portal.dists), method="single"))
	#plot(hclust(as.dist(portal.dists), method="average"))
	hc <- hclust(as.dist(portal.dists), method="centroid")
#	hc <- hclust(as.dist(portal.dists))
	# plot(hc)
	
	th <- 100
	
	repeat {
		
		portal.clusts <- list()
		height.lim <- TRUE
		for (i in 1:length(hc$height)) {
			if (height.lim && hc$height[i] < th) {
				portal.clust <- c()
				if (hc$merge[i,1] < 0) {
					portal.clust <- c(portal.clust, -(hc$merge[i,1]))
				} else {
					portal.clust <- c(portal.clust, portal.clusts[[hc$merge[i,1]]])
					portal.clusts[[hc$merge[i,1]]] <- NA
				}
				if (hc$merge[i,2] < 0) {
					portal.clust <- c(portal.clust, -(hc$merge[i,2]))
				} else {
					portal.clust <- c(portal.clust, portal.clusts[[hc$merge[i,2]]])
					portal.clusts[[hc$merge[i,2]]] <- NA
				}
				portal.clusts[[i]] <- portal.clust
			}
			else {
				height.lim <- FALSE
				if (hc$merge[i,1] < 0) {
					portal.clusts <- append(portal.clusts, list(- hc$merge[i,1]))
				}
				if (hc$merge[i,2] < 0) {
					portal.clusts <- append(portal.clusts, list(- hc$merge[i,2]))
				}
			}
		}
		
		portal.clusts <- portal.clusts[!is.na(portal.clusts)]
		
#		print(portal.clusts)

		# 1クラスタに7を超えるポータルが含まれる場合は、閾値を下げてクラスタを分割する
		if (th < 100) {
			if (max(unlist(lapply(portal.clusts, length))) > 7) {
				th <- th - 5
			} else {
				break
			} 
		}	
	
		# 総クラスタ数が８を超える場合は、閾値を下げてクラスタを統合する
		if (th > 100) {
			if (length(portal.clusts) > 8) {
				th <- th + 5
			} else {
				break
			} 
		}
		
		# 両方の閾値調整が必要な場合は、なにもしない
		if (th == 100) {
			if (max(unlist(lapply(portal.clusts, length))) > 6) {
				th <- th - 5
			} 
			if (length(portal.clusts) > 8) {
				th <- th + 5
			} 
			if (th == 100) {
				break
			}
		}
		
	}
	
	nclust <- length(portal.clusts)
	
	cat(nclust, "clusts", unlist(lapply(portal.clusts, length)), "\n")
	
	# portal.clusts
	
	
	# 各クラスタの重心を求める
	
	lnglat.center <- list()
	for (i in 1:nclust) {
		if (length(portal.clusts[[i]]) == 1) {
			lnglat.center[i] <- list(portals[mlcf.portals[portal.clusts[[i]]],c("lng", "lat")])
		} else {
			lnglat.center[i] <- list(apply(portals[mlcf.portals[portal.clusts[[i]]],c("lng", "lat")], 2, mean))
		}
	}
	
	lnglat.center <-  matrix(unlist(lnglat.center), byrow=T, ncol=2)
	colnames(lnglat.center) <- c("lng", "lat")
	
	# lnglat.center
	
	# クラスタの重心間の距離
	
	clust.dists <- c()
	for (i in 1:nclust) {
		clust.dists <- c(clust.dists, lnglat.dist(lnglat.center[i, "lng"], lnglat.center[i, "lat"], lnglat.center[, "lng"], lnglat.center[, "lat"]))
	}
	clust.dists <- matrix(clust.dists, byrow=T, ncol = nclust)
	
	# clust.dists
	
	# 頂点を含むクラスタを探す
	
	vertex.clusts <- c()
	for (v in c(1,2,3)) {
		for (i in 1:nclust) {
			if(is.element(v, portal.clusts[[i]])) {
				vertex.clusts[v] <- i
				break
			}
		}
	}
	
	# vertex.clusts
	
	best.patterns <- list()
	path.dists <- c()
	
	clust.pattern <- matrix(
		c(1,2,3, 2,3,1, 3,1,2), byrow=T, ncol=3)
	
	for (clust.pattern.index in 1:nrow(clust.pattern)) {
	
		start.clust <- vertex.clusts[clust.pattern[clust.pattern.index, 1]] 
		end.clust <- vertex.clusts[clust.pattern[clust.pattern.index, 3]] 
		message("start.clust ", start.clust, " end.clust ", end.clust)
		
		# クラスタの順番の一覧
		
		clust.path.pattern <- function(selected, rest) {
			if (length(rest) == 1) {
				return(c(selected, rest)) 
			}
			else {
				pattern <- c()
				for (i in 1:length(rest)) {
					pattern <- c(pattern, clust.path.pattern(c(selected, rest[i]), rest[-i]))
				}
				return(pattern)
			}
		}
		
		way.clusts <- c(1:nclust)[- c(start.clust, end.clust)]
		clust.path.patterns <- clust.path.pattern(c(), way.clusts)
		clust.path.patterns <- cbind(start.clust, matrix(clust.path.patterns, byrow=T, ncol=nclust - 2), end.clust)
		
	#	clust.path.patterns
	
		message(nrow(clust.path.patterns), "clust patterns")
		
		# クラスタの順番の総距離
		
		pattern.dists <- rep(0, nrow(clust.path.patterns))
		
		foreach (i = 1:nrow(clust.path.patterns)) %do% {
			pattern <- clust.path.patterns[i,]
			dist <- 0
			for(j in 2:nclust) {
				dist <- dist + clust.dists[pattern[j-1], pattern[j]] 
			}
			pattern.dists[i] <- dist
		}
		
		#	pattern.dists
		
		# クラスタの順番の総距離の一番短いもの
		
		min.pattern <- clust.path.patterns[which(pattern.dists == min(pattern.dists)),]
		
		min.pattern.clusts <- portal.clusts[min.pattern]
		
	#	for (i in 1:length(min.pattern.clusts)) {
	#		for (j in min.pattern.clusts[[i]]) {
	#			cat(portal.names[mlcf.portals[j]], "\n")
	#		}
	#		cat("\n")
	#	}
		
		# 順列を求める関数
		
		permu <- function(x, n) {
			iter <- function(n, x, ls) {
				if(n == 0){ ls }
				else{ sapply(x, function(xs) iter(n-1, x[x!=xs], rbind(ls, xs))) }
			}
			r <- iter(n, x, NULL)
			dim(r) <- c(n, length(r) / n)
			r
		}
		
		# クラスタ内のポータル順序
		# 先頭・最後が固定で、最短の中間の順序のものを選んでおく
		
		clust.portal.path.pattern <- list()
		for (clust.index in 1:nclust) {
			clust <- min.pattern.clusts[[clust.index]]
			path.patterns <- NA
			if (length(clust) == 1) {
				path.patterns <- matrix(clust, ncol=1)
			} else if (length(clust) == 2) {
				path.patterns <- matrix(c(clust[1], clust[2], clust[2], clust[1]), ncol=2)
				if (clust.index == 1) {
					path.patterns <- matrix(path.patterns[path.patterns[, 1] == clust.pattern[clust.pattern.index, 1], ], ncol=2)
				}
				if (clust.index == nclust) {
					path.patterns <- matrix(path.patterns[path.patterns[, 2] == clust.pattern[clust.pattern.index, 3], ], ncol=2)
				}
			} 
			else {
				pairs <- permu(clust, 2)

				if (clust.index == 1) {
					pairs <- pairs[, pairs[1, ] == clust.pattern[clust.pattern.index, 1]]
				}
				if (clust.index == nclust) {
					pairs <- pairs[, pairs[2, ] == clust.pattern[clust.pattern.index, 3]]
				}
				
				path.patterns <- c()
				for (pair.index in 1:ncol(pairs)) {
					pair <- pairs[, pair.index]
					rest <- setdiff(clust, pair)
					patterns <- permu(rest, length(rest))
					patterns <- rbind(pair[1], patterns, pair[2])
					n <- nrow(patterns)
					clust.portal.dists <- c()
					for (pattern.index in 1:ncol(patterns)) {
						pattern <- patterns[, pattern.index]
						dist <- 0
						for (portal.index in 2:n) {
							dist <- dist + portal.dists[pattern[portal.index-1], pattern[portal.index]]
						}
						clust.portal.dists <- c(clust.portal.dists, dist)
					}
					min.pattern <- patterns[, which(clust.portal.dists == min(clust.portal.dists))]
					path.patterns <- rbind(path.patterns, min.pattern)
				}
			}
			clust.portal.path.pattern[[clust.index]] <- path.patterns
		}
		
	#	clust.portal.path.pattern
		
		message(length(clust.portal.path.pattern), "clust portal path patterns")

		# ポータルの順番の一覧
	
		portal.path.pattern <- function(path.patterns) {
			if (length(path.patterns) == 1) {
				return(path.patterns[[1]])
			}
			
			self.patterns <- path.patterns[[1]]
			path.patterns[[1]] <- NULL
			child.pattern <- portal.path.pattern(path.patterns)
			
			npattern <- nrow(child.pattern)
			
			patterns <- c()
			for(i in 1:nrow(self.patterns)) {
				self.pattern <- self.patterns[i, ]
				pattern <- rep(self.pattern, npattern)
				dim(pattern) <- c(length(self.pattern), npattern)
				patterns <- rbind(patterns, cbind(t(pattern), child.pattern))
			}
			return (patterns)
		}
		
		portal.path.patterns <- portal.path.pattern(clust.portal.path.pattern)
		
		# portal.path.patterns
		
		# ポータルの順番の総距離の一番短いもの
	
		message(nrow(portal.path.patterns), "portal patterns")
	
		pattern.dists <- rep(0, nrow(portal.path.patterns))

		for (i in 1:nrow(portal.path.patterns)) {
			pattern <- portal.path.patterns[i,]
			dist <- 0
			for(j in 2:16) {
				dist <- dist + portal.dists[pattern[j-1], pattern[j]] 
			}
			pattern.dists[i] <- dist
		}
		
	#	pattern.dists
		
		best.patterns[[clust.pattern.index]] <- portal.path.patterns[which(pattern.dists == min(pattern.dists)),]
		path.dists[clust.pattern.index] <- min(pattern.dists)
	}
	
	
	# 経路パターン
	
	path.patterns <- list()
	
	for (pattern.index in 1:6) {
		
		best.pattern <- best.patterns[[(pattern.index + 1) %/% 2]] 
		if (pattern.index %% 2 == 0) {
			best.pattern <- rev(best.pattern)
		}
		
		path <- data.frame(posid = mlcf.posid[best.pattern], 
											 name = portal.names[mlcf.portals[best.pattern]],
											 stringsAsFactors = F)
	
		path.patterns[[pattern.index]] <- path
	}
}


# リンク先

link.def = c(
	1,2, 
	2,3,
	3,1,
	
	9,1,
	9,2,
	9,3,
	
	91, 1,
	91, 2,
	91, 9,
	92, 2,
	92, 3,
	92, 9,
	93, 3,
	93, 1,
	93, 9,
	
	911, 1,
	911, 2,
	911, 91,
	912, 2,
	912, 9,
	912, 91,
	913, 9,
	913, 1,
	913, 91,
	
	921, 2,
	921, 3,
	921, 92,
	922, 3,
	922, 9,
	922, 92,
	923, 9,
	923, 2,
	923, 92,
	
	931, 3,
	931, 1,
	931, 93,
	932, 1,
	932, 9,
	932, 93,
	933, 9,
	933, 3,
	933, 93
)

# リンク時の順序制約

constraints <- list(
	list(vertex=c(1,2,9),	inner=c(91,911,912,913)),
	list(vertex=c(2,3,9),	inner=c(92,921,922,923)),
	list(vertex=c(3,1,9),	inner=c(93,931,932,933)),
	list(vertex=c(1,2,91),	inner=c(911)),
	list(vertex=c(2,9,91),	inner=c(912)),
	list(vertex=c(1,9,91),	inner=c(913)),
	list(vertex=c(2,3,92),	inner=c(921)),
	list(vertex=c(3,9,92),	inner=c(922)),
	list(vertex=c(2,9,92),	inner=c(923)),
	list(vertex=c(3,1,93),	inner=c(931)),
	list(vertex=c(1,9,93),	inner=c(932)),
	list(vertex=c(3,9,93),	inner=c(933))
)

# リンク定義を双方向にする

link2 <- c()
for (i in seq(1,length(link.def)-1, 2)) {
	link2 <- c(link2, link.def[i], link.def[i+1], link.def[i+1], link.def[i])
} 
link <- matrix(link2, byrow=T, ncol=2)
colnames(link) <- c("p1", "p2")


check.path <- function(path) {

	# ポータル順制約を確認
	
	nposid <- length(path$posid)
	ng.portals <- c()
	
	for (constraint.index in 1:length(constraints)) {
		constraint.vertex <- constraints[[constraint.index]]$vertex
		constraint.inner <- constraints[[constraint.index]]$inner
		
		# 頂点部分の確保が何番目かを求める
		last.vertex <- max(which(is.element(path$posid, constraint.vertex)))
		
		# 最後の頂点確保よりも前に内部ポータルを確保する必要がある
		if (last.vertex < nposid) {
			check.inner <- is.element(constraint.inner, path$posid[(last.vertex+1):nposid])
			if (sum(check.inner) > 0) {
				ng <- constraint.inner[check.inner]
				ng.portals <- c(ng.portals, ng)
			}
		}
	}
	return(unique(ng.portals))
}

print.path <- function(path, ng.portals) {
	for (i in 1:nrow(path)) {
		cat(sprintf("%d\t%s\t%s\n",
								path$posid[i],
								path$name[i],
								ifelse(is.element(path$posid[i], ng.portals), "<NG>", "")))
	}
}

print.link.path <- function(path, ng.portals) {
	
	# リンク順を決める
	
	plist <- c()
	proc <- c()
	for (p1 in path$posid) {
		
		# すでに確保したポータルに対してリンクを張る
		
		if (is.null(plist)) {
			# 最初のポータルはリンクなし
			plist <- p1
			proc <- c(p1,NA)
		} else {
			# 確保したポータルから張れるリンク
			plist <- c(plist, p1)
			link.cand <- link[link[, "p1"] == p1, "p2"]
			p2s <- intersect(link.cand, plist)
			if (is.null(p2s)) {
				proc <- c(proc, p1,NA)
			} else {
				for (p2 in p2s) {
					proc <- c(proc, p1,p2)
				}
			}
		}
	}
	proc <- matrix(proc, byrow=T, ncol=2)
	colnames(proc) <- c("p1", "p2")
	
	# ポータル順（＋前ポータルとの距離）＋リンク順を出力
	
	last <- 0
	for (p1 in plist) {
		# 前ポータルとの距離を計算
		dist <- 0
		pos <- which(mlcf.posid == p1)
		if (last) {
			dist <- portal.dists[last, pos]
		}
		last <- pos
		
		# 順序に問題はないか
		is.ng <- ""
		if (is.element(p1, ng.portals)) {
			is.ng <- "<NG>\t"
		}
		cat(sprintf("%s%s\t%d\n", is.ng, path[path$posid == p1, "name"], round(dist)))
		
		for (p2 in proc[proc[, "p1"]==p1, "p2"]) {
			if (!is.na(p2)) {
				cat("\t", path[path$posid == p2, "name"], "\n")
			}
		}
	}
	
	cat("\n")
	
	# 必要な鍵数を出力
	
	for (p1 in mlcf.posid) {
		n <- sum(proc[, "p2"] == p1, na.rm=T)
		cat(path[path$posid == p1, "name"], "\t", n, "\n")
	}
	cat("\n")
	
	cat("[")
	cat(paste(paste(pls.min, collapse=","), collapse=""))
	cat(",{\"type\":\"polyline\",\"latLngs\":[")	
	
	initial <- T
	for (p1 in plist) {
		if (!initial) {
			cat(",")
		}
		initial <- F
		cat(sprintf("{\"lat\":%f,\"lng\":%f}", portals[mlcf.portals[p1 == mlcf.posid], "lat"], portals[mlcf.portals[p1 == mlcf.posid], "lng"]))
	}
	cat("],\"color\":\"#a24ac3\"}]\n")	
	cat("\n")	
}

polyline <- function(p1, p2) {
	return(sprintf('{"type":"polyline","latLngs":[{"lat":%f,"lng":%f},{"lat":%f,"lng":%f}],"color":"#ff0000"}', 
								 portals[p1,"lat"], portals[p1,"lng"], portals[p2,"lat"], portals[p2,"lng"]))
}

pls <- c(
	polyline(mlcf$vertexes[1], mlcf$vertexes[2]),
	polyline(mlcf$vertexes[2], mlcf$vertexes[3]),
	polyline(mlcf$vertexes[3], mlcf$vertexes[1]),
	polyline(inner$center, mlcf$vertexes[1]),
	polyline(inner$center, mlcf$vertexes[2]),
	polyline(inner$center, mlcf$vertexes[3]),
	polyline(inner$list1[[1]]$center, inner$center),
	polyline(inner$list1[[1]]$center, mlcf$vertexes[1]),
	polyline(inner$list1[[1]]$center, mlcf$vertexes[2]),
	polyline(inner$list2[[1]]$center, inner$center),
	polyline(inner$list2[[1]]$center, mlcf$vertexes[2]),
	polyline(inner$list2[[1]]$center, mlcf$vertexes[3]),
	polyline(inner$list3[[1]]$center, inner$center),
	polyline(inner$list3[[1]]$center, mlcf$vertexes[3]),
	polyline(inner$list3[[1]]$center, mlcf$vertexes[1]),
	
	polyline(inner$list1[[1]]$list1[1], inner$list1[[1]]$center),
	polyline(inner$list1[[1]]$list1[1], mlcf$vertexes[1]),
	polyline(inner$list1[[1]]$list1[1], mlcf$vertexes[2]),
	polyline(inner$list1[[1]]$list2[1], inner$list1[[1]]$center),
	polyline(inner$list1[[1]]$list2[1], mlcf$vertexes[2]),
	polyline(inner$list1[[1]]$list2[1], inner$center),
	polyline(inner$list1[[1]]$list3[1], inner$list1[[1]]$center),
	polyline(inner$list1[[1]]$list3[1], mlcf$vertexes[1]),
	polyline(inner$list1[[1]]$list3[1], inner$center),
	
	polyline(inner$list2[[1]]$list1[1], inner$list2[[1]]$center),
	polyline(inner$list2[[1]]$list1[1], mlcf$vertexes[2]),
	polyline(inner$list2[[1]]$list1[1], mlcf$vertexes[3]),
	polyline(inner$list2[[1]]$list2[1], inner$list2[[1]]$center),
	polyline(inner$list2[[1]]$list2[1], mlcf$vertexes[3]),
	polyline(inner$list2[[1]]$list2[1], inner$center),
	polyline(inner$list2[[1]]$list3[1], inner$list2[[1]]$center),
	polyline(inner$list2[[1]]$list3[1], mlcf$vertexes[2]),
	polyline(inner$list2[[1]]$list3[1], inner$center),
	
	polyline(inner$list3[[1]]$list1[1], inner$list3[[1]]$center),
	polyline(inner$list3[[1]]$list1[1], mlcf$vertexes[3]),
	polyline(inner$list3[[1]]$list1[1], mlcf$vertexes[1]),
	polyline(inner$list3[[1]]$list2[1], inner$list3[[1]]$center),
	polyline(inner$list3[[1]]$list2[1], mlcf$vertexes[1]),
	polyline(inner$list3[[1]]$list2[1], inner$center),
	polyline(inner$list3[[1]]$list3[1], inner$list3[[1]]$center),
	polyline(inner$list3[[1]]$list3[1], mlcf$vertexes[3]),
	polyline(inner$list3[[1]]$list3[1], inner$center)
)

pls.min <- c(
	polyline(mlcf$vertexes[1], mlcf$vertexes[2]),
	polyline(mlcf$vertexes[2], mlcf$vertexes[3]),
	polyline(mlcf$vertexes[3], mlcf$vertexes[1]),
	polyline(inner$center, mlcf$vertexes[1]),
	polyline(inner$center, mlcf$vertexes[2]),
	polyline(inner$center, mlcf$vertexes[3]))




# パスパターン（ポータル順とエラーポータル）

ng.portals <- list()

for (pattern.index in 1:6) {

	ng <- check.path(path.patterns[[pattern.index]])
	ng.portals[[pattern.index]] <- ifelse(is.null(ng), NA, ng)

}

if (!debug) {
	sink(proc.filename)
}

# ポータル一覧

for (i in 1:length(mlcf.posid)) {
	cat(sprintf("%d\t%s\n", mlcf.posid[i], portal.names[mlcf.portals[i]]))
}
cat("\n")	


vertex.pattern <- matrix(
	c(1,2,3, 3,2,1, 2,3,1, 1,3,2, 3,1,2, 2,1,3), byrow=T, ncol=3)

for (pattern.index in 1:6) {

	plist <- path.patterns[[pattern.index]]$posid

	dist <- 0
	last <- 0
	for (p1 in plist) {
		# 前ポータルとの距離を計算
		pos <- which(mlcf.posid == p1)
		if (last) {
			dist <- dist + portal.dists[last, pos]
		}
		last <- pos
	}
	
	cat(sprintf("%s | %s | %s (%d) %s\n",
							portal.names[mlcf.portals[vertex.pattern[pattern.index, 1]]],
							portal.names[mlcf.portals[vertex.pattern[pattern.index, 2]]],
							portal.names[mlcf.portals[vertex.pattern[pattern.index, 3]]],
							round(dist),
							ifelse(is.na(ng.portals[[pattern.index]]), "", "<NG>")))
	
}
cat("\n")	

# リンク一覧

cat(paste("[", paste(pls, collapse=","), "]", collapse=""))
cat("\n")
cat("\n\n")


for (pattern.index in 1:6) {
	print.path(path.patterns[[pattern.index]], ng.portals[[pattern.index]])
	cat("\n")
}
cat("\n")
cat("\n")

for (pattern.index in 1:6) {
	print.link.path(path.patterns[[pattern.index]], ng.portals[[pattern.index]])
}
cat("\n")

if (!debug) {
	sink()
}

