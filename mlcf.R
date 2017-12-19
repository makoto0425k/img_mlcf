source("parallel.R")
source("progress.R")

args <- commandArgs(trailingOnly = T)

cpu <- 4
chunk <- 4
layer <- 4
limit <- c(0,0,0,0,0)
far <- 1000
clear <- F
main <- c()

args <- commandArgs(trailingOnly=TRUE)
if (length(args)) {
	for (i in 1:length(args)) {
		if (substring(args[i], 1, 4) == "-cpu") {
			cpu <- chunk <- as.numeric(substring(args[i], 5))
		} else if (args[i] == "-4") {
			layer <- 4
			far <- 1000
			limit <- c(0,0,0,0,0)
		} else if (args[i] == "-5") {
			layer <- 5
			far <- 1500
			limit <- c(0,4,4,4,4)
		} else if (substring(args[i], 1, 2) == "-4" && nchar(args[i]) == 7) {
			layer <- 4
			far <- 1000
			limit <- as.numeric(substring(args[i], 3:7, 3:7))
		} else if (substring(args[i], 1, 2) == "-5" && nchar(args[i]) == 7) {
			layer <- 5
			far <- 1500
			limit <- as.numeric(substring(args[i], 3:7, 3:7))
		} else if (substring(args[i], 1, 2) == "-f") {
			far <- as.numeric(substring(args[i], 3))
		} else if (substring(args[i], 1, 2) == "-n") {
			chunk <- as.numeric(substring(args[i], 3))
		} else if (substring(args[i], 1, 5) == "-main") {
			main <- 1:as.numeric(substring(args[i], 6))
		} else if (args[i] == "-clear") {
			clear <- TRUE
		} else {
			filename <- args[i]
		}
	}
}

if (length(args) == 0) {
	source("defaultargs.R")
}

message("file ", filename, "  main ", ifelse(is.null(main), 0, max(main)))
message("layer ", layer, "  far ", far)
message("cpu ", cpu, "  chunk ", chunk)


# 完全N重の内側のポータル数
limit.0 <- c(0,1,4,13,40)
# 均一N重のための内側ポータル数の上限、limit分まで不使用ポータルを含んでOK
limit.1 <- limit.0 + limit 


# 見つかった多重パターンをgooglemapに表示する
# テンプレートを読み込んでtextareaあたりに結果を埋め込む
# テンプレートには、html表示時にtextareaの内容をパースしてgooglemapに情報を表示するようなjavascriptを埋め込んである

mlcf.html <- scan("mlcf.html", what = character(), sep = "\n", blank.lines.skip = F)
portal.info.line <- "          <textarea id=\"triangle_info\"></textarea>"
portal.info.line.no <- which(mlcf.html == portal.info.line)
range.info.line <- "          <input type=\"hidden\" id=\"map_range\" value=\"\">"	
range.info.line.no <- which(mlcf.html == range.info.line)

print.map <- function(mlcfs, layer) {
	
	info <- sapply (mlcfs, function(mlcf) {
		vertex1 <- mlcf$vertexes[1]
		vertex2 <- mlcf$vertexes[2]
		vertex3 <- mlcf$vertexes[3]
		center <- mlcf$inner[[1]]$center
		
		paste(
			sprintf("(%d) %s | %s | %s\t", mlcf$count, portal.list[vertex1, "name"], portal.list[vertex2, "name"], portal.list[vertex3, "name"], portals[center,"lat"], portals[center,"lng"]),
			sprintf("%f,%f\t", portals[center,"lat"], portals[center,"lng"]),
			sprintf("%f,%f,%f,%f\t", portals[vertex1,"lat"], portals[vertex1,"lng"], portals[vertex2,"lat"], portals[vertex2,"lng"]),
			sprintf("%f,%f,%f,%f\t", portals[vertex2,"lat"], portals[vertex2,"lng"], portals[vertex3,"lat"], portals[vertex3,"lng"]),
			sprintf("%f,%f,%f,%f\t", portals[vertex3,"lat"], portals[vertex3,"lng"], portals[vertex1,"lat"], portals[vertex1,"lng"]), 
			sprintf("%f,%f,%f\n", 
							lnglat.dist(portals[vertex1,"lng"], portals[vertex1,"lat"], portals[vertex2,"lng"], portals[vertex2,"lat"] ),
							lnglat.dist(portals[vertex2,"lng"], portals[vertex2,"lat"], portals[vertex3,"lng"], portals[vertex3,"lat"] ),
							lnglat.dist(portals[vertex3,"lng"], portals[vertex3,"lat"], portals[vertex1,"lng"], portals[vertex1,"lat"] )))
	})

	range.info <- sprintf("%f,%f,%f,%f", latlng.min[1], latlng.min[2], latlng.max[1], latlng.max[2])

	mlcf.html[portal.info.line.no] <- paste("          <textarea id=\"triangle_info\">", paste(1:length(info), info, collapse="") , "</textarea>", sep="")
	mlcf.html[range.info.line.no ] <- paste("          <input type=\"hidden\" id=\"map_range\" value=\"", range.info , "\">", sep="")
	
	if (layer == 4) {	keyword <- "mlcf4" }
	if (layer == 5) {	keyword <- "mlcf5" }
	
	out <- file(paste(filename, ".", keyword, ".map.html", sep=""), "w", encoding="utf-8")
	writeLines(mlcf.html, out, sep="\n")
	close(out)
	
	invisible(0)
}


# 完全（均一）多重を構成できるかどうかをチェック
# 内側ポータルを１個選び、頂点＋選んだポータルで構成される３つのフィールドが全て均一多重（多重度が１低い）かを調べる

check.mlcf <- function(layer, sele, c1, c2, c3, limit.0, limit.1) {
	
	layer <- layer - 1	  
	
	r <- list();

	# 頂点と他ポータルに対する方向をまとめて計算
	vertex <- portals[c1,]
	r1 <- (
			matrix(c((portals[sele, "lng"] - vertex["lng"]), (portals[sele, "lat"] - vertex["lat"])), byrow=F, nrow=sum(sele)) 
				%*% 
			matrix(c((portals[sele, "lat"] - vertex["lat"]), (vertex["lng"] - portals[sele, "lng"])), byrow=T, ncol=sum(sele)) 
		) < 0
	
	vertex <- portals[c2,]
	r2 <- (
			matrix(c((portals[sele, "lng"] - vertex["lng"]), (portals[sele, "lat"] - vertex["lat"])), byrow=F, nrow=sum(sele)) 
				%*% 
			matrix(c((portals[sele, "lat"] - vertex["lat"]), (vertex["lng"] - portals[sele, "lng"])), byrow=T, ncol=sum(sele)) 
		) < 0
	
	vertex <- portals[c3,]
	r3 <- (
			matrix(c((portals[sele, "lng"] - vertex["lng"]), (portals[sele, "lat"] - vertex["lat"])), byrow=F, nrow=sum(sele)) 
				%*% 
			matrix(c((portals[sele, "lat"] - vertex["lat"]), (vertex["lng"] - portals[sele, "lng"])), byrow=T, ncol=sum(sele)) 
		) < 0
	
	sele.c1 <- which(sele) == c1
	sele.c2 <- which(sele) == c2
	sele.c3 <- which(sele) == c3

	# 対象ポータルを中心にできるかをチェック
	for (c4 in which(sele)) {

		sele.c4 <- which(sele) == c4

		sele12in <- r1[, sele.c4] == r1[sele.c2, sele.c4] & r2[, sele.c4] == r2[sele.c1, sele.c4]
		sele12in[sele.c1 | sele.c2 | sele.c4] <- T
		if ((sum(sele12in) - 3) < limit.0[layer] || (sum(sele12in) - 3) > limit.1[layer]) { next }

		sele23in <- r2[, sele.c4] == r2[sele.c3, sele.c4] & r3[, sele.c4] == r3[sele.c2, sele.c4]
		sele23in[sele.c2 | sele.c3 | sele.c4] <- T
		if ((sum(sele23in) - 3) < limit.0[layer] || (sum(sele23in) - 3) > limit.1[layer]) { next }
		
		sele31in <- r3[, sele.c4] == r3[sele.c1, sele.c4] & r1[, sele.c4] == r1[sele.c3, sele.c4]
		sele31in[sele.c3 | sele.c1 | sele.c4] <- T
		if ((sum(sele31in) - 3) < limit.0[layer] || (sum(sele31in) - 3) > limit.1[layer]) { next }

		if (layer == 2) {
			sele12in[sele.c1 | sele.c2 | sele.c4] <- F
			sele23in[sele.c2 | sele.c3 | sele.c4] <- F
			sele31in[sele.c3 | sele.c1 | sele.c4] <- F

			sele12 <- sele
			sele12[sele] <- sele12in
			sele23 <- sele
			sele23[sele] <- sele23in
			sele31 <- sele
			sele31[sele] <- sele31in
			
			r <- append(r, list(list(center=c4, list1=which(sele12), list2=which(sele23), list3=which(sele31))))
		}	else {
			sele12 <- sele
			sele12[sele] <- sele12in
			sele23 <- sele
			sele23[sele] <- sele23in
			sele31 <- sele
			sele31[sele] <- sele31in
			
			p1 <- check.mlcf(layer, sele12, c1, c2, c4, limit.0, limit.1)
			if (is.null(p1)) { next; }
			p2 <- check.mlcf(layer, sele23, c2, c3, c4, limit.0, limit.1)
			if (is.null(p2)) { next; }
			p3 <- check.mlcf(layer, sele31, c3, c1, c4, limit.0, limit.1)
			if (is.null(p3)) { next; }
			
			r <- append(r, list(list(center=c4, list1=p1, list2=p2, list3=p3)))
		}
	}
	
	if (length(r) == 0) {
		return(c())
	}
	
	return(r)
}


# ポータル一覧の読み込み

portal.list <- read.csv(paste(filename, ".txt", sep=""), header=F, col.names=c("name","lat","lng"), fileEncoding="UTF-8", encoding="UTF-8", stringsAsFactors = F)
portals <- as.matrix(portal.list[, c("lat","lng")])

latlng.min <- apply(portals, 2, min)
latlng.max <- apply(portals, 2, max)


# 距離マトリクスを作成

source("lnglatdist.R")

portal.count <- nrow(portals)

parallel.start(cpu)

dists <- foreach (i = 1:portal.count, .combine = rbind, .inorder=T) %dopar% {
	lnglat.dist(portals[i, "lng"], portals[i, "lat"], portals[, "lng"], portals[, "lat"])
}

parallel.end()


# 残り処理時間を計算するために、全体の計算量を求めておく

if (length(main) == 0) {
	main <- 1:portal.count
}
main.portal.count <- length(main)

cf.count <- foreach (c1 = main[1:(length(main) - 2)], .combine='c') %do% {
	
	portals.target.2 <- dists[c1, ] < far
	portals.target.2[1:c1] <- FALSE
	
	if (sum(portals.target.2) == 0) {
		return(0)
	}
	
	return (sum(sapply (which(portals.target.2[1:(portal.count - 1)]), function(c2) {
		
		portals.target.3 <- (dists[c2, ] < far) & portals.target.2
		portals.target.3[1:c2] <- FALSE
		
		return(sum(portals.target.3))
	})))
	
}


# 実行経過があればそれを読み込み
# （途中まで実行したところで、何らかの理由でプログラム終了したときに、その時点から処理を再開できるようにしてある）

mlcfs <- list()
c1.proc <- 0

if (!clear) {
	if (layer == 4 && file.exists(paste(filename, ".mlcf4.rdat", sep=""))) {
		load(paste(filename, ".mlcf4.rdat", sep=""))
	}
	if (layer == 5 && file.exists(paste(filename, ".mlcf5.rdat", sep=""))) {
		load(paste(filename, ".mlcf5.rdat", sep=""))
	}
	
	if (file.exists(paste(filename, ".proc.rdat", sep=""))) {
		load(paste(filename, ".proc.rdat", sep=""))
	}
}


# 残り時間（推定）・経過時間の表示のための初期化

calc.vol <- sum(cf.count)
skip.vol <- 0
if (c1.proc > 0) {
	skip.vol <- sum(cf.count[1:c1.proc])
}
elapsed.vol <- skip.vol
progress.start(calc.vol)


# ここからが主処理

parallel.start(cpu)
parallel.export <- c("check.mlcf")

while(c1.proc < main.portal.count) {

  # 経過状況を表示
	comment <- sprintf("%s %d %d/%d mlcf:%d", filename, layer, c1.proc, main.portal.count, length(mlcfs))
	calc.vol <- 0
	if (c1.proc > 0) {
		calc.vol <- sum(cf.count[1:c1.proc])
	}
	progress.report(calc.vol, skip = skip.vol, comment=comment)
	
	c1.chunk <- (c1.proc + 1):min(c1.proc + chunk, main.portal.count)

	mlcfs1 <- foreach (c1 = c1.chunk, .combine = append, .export=parallel.export) %dopar% {

	  # １地点目選択
	  
	  c1.c2 <- dists[c1, ] < far
		c1.c2[1:c1] <- FALSE
		c1.c2[portal.count] <- FALSE

		mlcfs2 <- list()
		for (c2 in which(c1.c2)) {

		  # ２地点目選択
		  
		  vertex1 <- portals[c1, ]
			vertex2 <- portals[c2, ]
			
			# ２地点からの距離が近いポータルだけを選択
			# booleanベクトル、長さはポータル数
			sele12 <- dists[c1,] < far & dists[c2,] < far

			# 数が少なければ抜ける
			if ((sum(sele12) - 3) < limit.0[layer]) {	next	}
			
			# c1・c2 に対する方向
			# 1列N行、行数はsele12がtrueの（２地点からの距離が近いポータルの）個数
			r12 <- ( 
				matrix(c((portals[sele12, "lng"] - vertex1["lng"]), (portals[sele12, "lat"] - vertex1["lat"])), byrow=F, nrow=sum(sele12))
				%*% 
					matrix(c((vertex2["lat"] - vertex1["lat"]), (vertex1["lng"] - vertex2["lng"])), byrow=T, nrow=2) 
			) < 0
			
			# c1と他ポータルに対する方向をまとめて計算
			# N列N行、行数・列数はsele12がtrueの（２地点からの距離が近いポータルの）個数
			r1 <- (
				matrix(c((portals[sele12, "lng"] - vertex1["lng"]), (portals[sele12, "lat"] - vertex1["lat"])), byrow=F, nrow=sum(sele12)) 
				%*% 
					matrix(c((portals[sele12, "lat"] - vertex1["lat"]), (vertex1["lng"] - portals[sele12, "lng"])), byrow=T, ncol=sum(sele12)) 
			) < 0
			
			# c2と他ポータルに対する方向をまとめて計算
			# N列N行、行数・列数はsele12がtrueの（２地点からの距離が近いポータルの）個数
			r2 <- (
				matrix(c((portals[sele12, "lng"] - vertex2["lng"]), (portals[sele12, "lat"] - vertex2["lat"])), byrow=F, nrow=sum(sele12)) 
				%*% 
					matrix(c((portals[sele12, "lat"] - vertex2["lat"]), (vertex2["lng"] - portals[sele12, "lng"])), byrow=T, ncol=sum(sele12))
			) < 0
			
			# sele12の中でc1・c2は何番目かを判定するためのベクトル
			# booleanのベクトル、長さはsele12がtrueの（２地点からの距離が近いポータルの）個数
			sele.c1 <- which(sele12) == c1
			sele.c2 <- which(sele12) == c2
			
			# c3の候補		
			c1.c2.c3 <- sele12
			c1.c2.c3[1:c2] <- FALSE

			ra <- list()
			for (c3 in which(c1.c2.c3)) {
				
				# ３地点の内側
				sele.c3 <- which(sele12) == c3
				sele123 <- 
				  r12[, 1] == r12[sele.c3, 1] &  # 線分 c1・c2 で c3と同じ方向にあるもの
				  r1[, sele.c3] == r1[sele.c2, sele.c3] &  # 線分 c1・c3 で c2と同じ方向にあるもの
				  r2[, sele.c3] == r2[sele.c1, sele.c3]    # 線分 c2・c3 で c1と同じ方向にあるもの

				# sele123 は booleanベクトル、長さはsele12がtrueの（２地点からの距離が近いポータルの）個数
				
				sele123[sele.c1 | sele.c2 | sele.c3] <- T
				
				# 数が少なければ、または、数が多ければ、抜ける
				if ((sum(sele123) - 3) < limit.0[layer] || (sum(sele123) - 3) > limit.1[layer]) {	next	}
				
				# sele123を全体に拡張
				sele <- sele12
				sele[sele12] <- sele123
				# sele は booleanのベクトル、流さはポータル数、３点頂点および内部がtrue
				
				# 残ったポータル一覧で多重を構成できるか
				mlcf <- check.mlcf(layer, sele, c1, c2, c3, limit.0, limit.1)
				if (is.null(mlcf)) {
					next
				}
				
				ra <- append(ra, list(list(vertexes = c(c1,c2,c3), count = sum(sele) - 3, inner = mlcf)))
			}
			
			if (length(ra) == 0) {
				# return(NULL)
				next
			}
			
			# return(ra)
			mlcfs2 <- append(mlcfs2, ra)
		}
		
		return(mlcfs2)
	}

	parallel.export <- NULL
	
	if (length(mlcfs1)) {
		mlcfs1 <- mlcfs1[!sapply(mlcfs1, is.null)]
		if (length(mlcfs1) > 0) {
			mlcfs <- append(mlcfs, mlcfs1)
			if (layer == 4) {	ext <- ".mlcf4.rdat"	}
			if (layer == 5) {	ext <- ".mlcf5.rdat"	}
			save(mlcfs, file = paste(filename, ext, sep=""))
			print.map(mlcfs, layer)
		}
	}
	c1.proc <- c1.proc + length(c1.chunk)
	save(c1.proc, file = paste(filename, ".proc.rdat", sep=""))
}

comment <- sprintf("%s %d %d/%d mlcf:%d", filename, layer, c1.proc, main.portal.count, length(mlcfs))
progress.report(sum(cf.count), skip = skip.vol, comment=comment)
progress.end(skip = skip.vol)

parallel.end()

if (length(mlcfs)) {
	print.map(mlcfs, layer)
}

if (file.exists(paste(filename, ".proc.rdat", sep=""))) {
	file.remove(paste(filename, ".proc.rdat", sep=""))
}

invisible(0)
