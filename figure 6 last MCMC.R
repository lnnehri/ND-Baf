set.seed(1)

### 1) Parametreler (kolay ayarlanır)
grid_size <- 100
num_steps <- 2000
n_lys <- 300

# Bölgeler
nucleus_region <- c(40, 60, 40, 60)      # x1,x2,y1,y2
boundary_size <- 12                      # perifer tanımı (plazma membrana yakın)

# Koşullar: ND, ND+Baf gibi "durumlar" için parametre seti
params <- list(
  ND = list(
    p_alk_base = 0.15,     # alkalin lizozom olma baz olasılığı
    ca_gain = 0.04,        # exocytosis olunca Ca artışı
    p_exo_base = 0.002,    # baz exocytosis
    fao = 0.45             # hayatta kalma/enerji (0-1)
  ),
  ND_Baf = list(
    p_alk_base = 0.55,
    ca_gain = 0.09,
    p_exo_base = 0.010,
    fao = 0.60
  )
)

# Fenotip eşikleri
exo_threshold_for_elong <- 25   # total exocytosis eşiği
ca_threshold_for_elong  <- 1.5  # Ca birikim eşiği
death_prob_base <- 0.0005       # baz ölüm olasılığı (zaman adımı başına)

### 2) Bölge kontrol fonksiyonları
is_in_rect <- function(x, y, r) r[1] <= x && x < r[2] && r[3] <= y && y < r[4]

is_peripheral <- function(x, y, grid_size, boundary_size) {
  (x < boundary_size) || (x >= grid_size - boundary_size) ||
    (y < boundary_size) || (y >= grid_size - boundary_size)
}

### 3) Lizozomları başlat
# positions: 0-based (x,y)
pos <- cbind(
  sample.int(grid_size, n_lys, replace = TRUE) - 1,
  sample.int(grid_size, n_lys, replace = TRUE) - 1
)
colnames(pos) <- c("x","y")

# Lizozom tipi: 1=alkalin, 0=asidik (başlangıç)
alk <- rbinom(n_lys, 1, prob = 0.2)

### 4) Basit diffusion (nükleus bölgesinde kalma eğilimi gibi ekstra kurallar eklenebilir)
move_lysosomes <- function(pos, grid_size) {
  for (i in seq_len(nrow(pos))) {
    dir <- sample(c("up","down","left","right"), 1)
    x <- pos[i,1]; y <- pos[i,2]
    if (dir == "up" && x > 0) x <- x - 1
    if (dir == "down" && x < grid_size-1) x <- x + 1
    if (dir == "left" && y > 0) y <- y - 1
    if (dir == "right" && y < grid_size-1) y <- y + 1
    pos[i,] <- c(x,y)
  }
  pos
}

### 5) Tek koşul simülatörü (ND veya ND_Baf)
simulate_condition <- function(cond_name, params, pos0, alk0) {
  p <- params[[cond_name]]
  
  pos <- pos0
  alk <- alk0
  
  # hücre-düzeyi değişkenler
  Ca <- 0
  total_exo <- 0
  phenotype <- "rounded"   # rounded / elongated / dead
  
  # kayıt
  traj <- data.frame(step=integer(), Ca=double(), exo=integer(),
                     phen=character(), frac_periph_lys=double(),
                     frac_alk=double(), stringsAsFactors = FALSE)
  
  for (step in 0:(num_steps-1)) {
    
    if (phenotype == "dead") {
      traj <- rbind(traj, data.frame(step=step, Ca=Ca, exo=total_exo,
                                     phen=phenotype, frac_periph_lys=NA,
                                     frac_alk=mean(alk), stringsAsFactors=FALSE))
      next
    }
    
    # hareket
    pos <- move_lysosomes(pos, grid_size)
    
    # Koşula göre alkalinleşme (çok basit: her adımda küçük bir "flip" ihtimali)
    # ND+Baf: alk olma olasılığı daha yüksek
    flip_to_alk <- runif(n_lys) < (p$p_alk_base * 0.002)
    alk[flip_to_alk] <- 1
    
    # periferik lizozomlar
    periph <- mapply(is_peripheral, pos[,1], pos[,2],
                     MoreArgs=list(grid_size=grid_size, boundary_size=boundary_size))
    
    # exocytosis olasılığı: perifer + alkalin + baz
    # (istersen nucleus uzaklığı / LAMP1 protrusion vs eklenir)
    p_exo <- p$p_exo_base +
      0.03 * (periph & alk==1) +   # ana mekanizma
      0.005 * (periph & alk==0)
    
    # exocytosis olayları
    exo_events <- rbinom(n_lys, 1, prob = pmin(p_exo, 0.9))
    exo_count <- sum(exo_events)
    
    total_exo <- total_exo + exo_count
    Ca <- Ca + exo_count * p$ca_gain
    
    # FAO düşükse ölüm olasılığı artar (çok basit)
    p_death <- death_prob_base + (1 - p$fao) * 0.002
    if (runif(1) < p_death) phenotype <- "dead"
    
    # Elongation kararı (eşik tabanlı)
    if (phenotype != "dead") {
      if (total_exo >= exo_threshold_for_elong && Ca >= ca_threshold_for_elong) {
        phenotype <- "elongated"
      }
    }
    
    traj <- rbind(traj, data.frame(
      step = step,
      Ca = Ca,
      exo = total_exo,
      phen = phenotype,
      frac_periph_lys = mean(periph),
      frac_alk = mean(alk),
      stringsAsFactors = FALSE
    ))
    
    # elongated olduktan sonra istersen kır:
    # if (phenotype == "elongated") break
  }
  
  list(traj=traj, final=list(Ca=Ca, exo=total_exo, phen=phenotype))
}

### 6) Koşulları çalıştır
res_ND     <- simulate_condition("ND", params, pos, alk)
res_ND_Baf <- simulate_condition("ND_Baf", params, pos, alk)

### 7) Hızlı özet + basit plot
cat("ND final phenotype:", res_ND$final$phen, "| exo:", res_ND$final$exo, "| Ca:", round(res_ND$final$Ca,2), "\n")
cat("ND+Baf final phenotype:", res_ND_Baf$final$phen, "| exo:", res_ND_Baf$final$exo, "| Ca:", round(res_ND_Baf$final$Ca,2), "\n")

# Ca traj plot
plot(res_ND$traj$step, res_ND$traj$Ca, type="l", xlab="step", ylab="Ca (proxy)")
lines(res_ND_Baf$traj$step, res_ND_Baf$traj$Ca)

# phenotype zaman içinde (basit)
table(res_ND$traj$phen)
table(res_ND_Baf$traj$phen)
###coklu hucreyle
## Same model, more cells + a few tweaks to get multiple hotspots
set.seed(2)

## 1) Parameters
G <- 90                 # tissue grid size
N <- 220                # <-- more cells
Tmax <- 500

p_move_E <- 0.55
p_move_R <- 0.10

boundary_size <- 12

p_alk_base <- 0.01
boost_edge <- 0.06
boost_contact <- 0.10     # slightly higher -> more emission

p_emit_if_edge_alk_contact <- 0.45  # slightly higher -> more hotspots
signal_decay <- 0.93
signal_diffuse <- 0.18

E_threshold <- 1.0        # slightly lower -> easier R->E
D_prob_base <- 0.001

## 2) Helpers
clip <- function(x, lo=0, hi=G-1) pmax(lo, pmin(hi, x))
nbr4 <- function(x,y) rbind(c(x-1,y), c(x+1,y), c(x,y-1), c(x,y+1))

is_peripheral <- function(x, y, G, b) {
  (x < b) || (x >= G - b) || (y < b) || (y >= G - b)
}

## 3) Init agents
agents <- data.frame(
  id = 1:N,
  x = sample(0:(G-1), N, TRUE),
  y = sample(0:(G-1), N, TRUE),
  state = rep("R", N),
  stringsAsFactors = FALSE
)

H <- matrix(0, nrow=G, ncol=G)

## 4) Contact check
has_contact <- function(ax, ay, others){
  nb <- nbr4(ax, ay)
  for(k in 1:nrow(nb)){
    if(any(others$x == nb[k,1] & others$y == nb[k,2])) return(TRUE)
  }
  FALSE
}

## 5) Edge-alk draw (summary proxy)
draw_edge_alk <- function(in_contact){
  p <- p_alk_base + boost_edge + if(in_contact) boost_contact else 0
  runif(1) < pmin(p, 0.95)
}

## 6) Signal update
update_signal_field <- function(H){
  H <- H * signal_decay
  Hnew <- H
  for(i in 2:(G-1)){
    for(j in 2:(G-1)){
      lap <- (H[i-1,j] + H[i+1,j] + H[i,j-1] + H[i,j+1] - 4*H[i,j])
      Hnew[i,j] <- H[i,j] + signal_diffuse * lap
    }
  }
  Hnew
}

## 7) Movement
move_agent <- function(x,y,state,H){
  if(state=="D") return(c(x,y))
  p_move <- if(state=="E") p_move_E else p_move_R
  if(runif(1) > p_move) return(c(x,y))
  
  nb <- nbr4(x,y)
  nb[,1] <- clip(nb[,1]); nb[,2] <- clip(nb[,2])
  
  if(state=="E"){
    vals <- mapply(function(xx,yy) H[xx+1,yy+1], nb[,1], nb[,2])
    nb[which.max(vals), ]
  } else {
    nb[sample(1:4,1), ]
  }
}

## 8) Run
log_df <- data.frame(step=integer(), nR=integer(), nE=integer(), nD=integer())

for(t in 1:Tmax){
  
  alive <- agents[agents$state!="D", ]
  
  ## emissions
  for(i in 1:nrow(alive)){
    id <- alive$id[i]
    ax <- alive$x[i]; ay <- alive$y[i]
    others <- alive[alive$id!=id, ]
    
    in_contact <- has_contact(ax, ay, others)
    edge_alk <- draw_edge_alk(in_contact)
    
    # also allow peripheral bias at tissue boundary (optional)
    periph_tissue <- is_peripheral(ax, ay, G, boundary_size)
    
    if(edge_alk && in_contact && runif(1) < p_emit_if_edge_alk_contact){
      H[ax+1, ay+1] <- H[ax+1, ay+1] + 0.55
      if(periph_tissue) H[ax+1, ay+1] <- H[ax+1, ay+1] + 0.15
    }
  }
  
  H <- update_signal_field(H)
  
  ## transitions
  for(i in 1:N){
    if(agents$state[i]=="D") next
    ax <- agents$x[i]; ay <- agents$y[i]
    localH <- H[ax+1, ay+1]
    
    if(runif(1) < D_prob_base){
      agents$state[i] <- "D"; next
    }
    if(agents$state[i]=="R" && localH >= E_threshold) agents$state[i] <- "E"
  }
  
  ## movement with collision avoidance
  for(i in 1:N){
    if(agents$state[i]=="D") next
    newxy <- move_agent(agents$x[i], agents$y[i], agents$state[i], H)
    
    occ <- agents$state!="D"
    if(!any(agents$x[occ]==newxy[1] & agents$y[occ]==newxy[2])){
      agents$x[i] <- newxy[1]; agents$y[i] <- newxy[2]
    }
  }
  
  log_df <- rbind(log_df, data.frame(
    step=t,
    nR=sum(agents$state=="R"),
    nE=sum(agents$state=="E"),
    nD=sum(agents$state=="D")
  ))
}

## 9) Plots
plot(log_df$step, log_df$nE, type="l", xlab="step", ylab="# elongated (E)")

image(t(H[nrow(H):1, ]), axes=FALSE, main="Final heterogeneity signal field (more cells)")
box()


#1 hucre
## 0) MINIMAL grid-based ABM (single-cell + population) — sequential script
set.seed(1)

## 1) Parameters
G <- 80                 # tissue grid size (GxG)
N <- 60                 # number of cells (agents)
Tmax <- 400             # steps

p_move_E <- 0.55        # move prob if elongated
p_move_R <- 0.10        # move prob if rounded

# "edge" inside a cell (how many tiles from cell center count as "edge")
cell_radius <- 3
edge_band <- 1

# alkalinization baseline + boosts
p_alk_base <- 0.01
boost_edge <- 0.06
boost_contact <- 0.08

# heterogeneity signal creation
p_emit_if_edge_alk_contact <- 0.35
signal_decay <- 0.92
signal_diffuse <- 0.15

# phenotype transitions
E_threshold <- 1.2      # if local signal exceeds -> R->E
D_prob_base <- 0.001

## 2) Helpers: bounds + neighborhood
clip <- function(x, lo=0, hi=G-1) pmax(lo, pmin(hi, x))

nbr4 <- function(x,y){
  rbind(c(x-1,y), c(x+1,y), c(x,y-1), c(x,y+1))
}

## 3) Initialize population (agent positions, states)
agents <- data.frame(
  id = 1:N,
  x = sample(0:(G-1), N, TRUE),
  y = sample(0:(G-1), N, TRUE),
  state = rep("R", N),          # R, E, D
  stringsAsFactors = FALSE
)

# tissue-level signal field H(x,y)
H <- matrix(0, nrow=G, ncol=G)

## 4) Contact map (cell-cell contact if any other cell in 4-neighborhood)
has_contact <- function(ax, ay, agents_alive){
  nb <- nbr4(ax, ay)
  for(k in 1:nrow(nb)){
    if(any(agents_alive$x == nb[k,1] & agents_alive$y == nb[k,2])) return(TRUE)
  }
  FALSE
}

## 5) Single-cell internal "mini-grid" (lysosome edge alkalinization proxy)
# We don't simulate each lysosome; we simulate a probability summary:
# edge_alk = TRUE/FALSE each time step per cell, driven by edge bias + contact bias.
draw_edge_alk <- function(in_contact){
  p <- p_alk_base + boost_edge + if(in_contact) boost_contact else 0
  runif(1) < pmin(p, 0.95)
}

## 6) Signal diffusion + decay on tissue grid
update_signal_field <- function(H){
  H <- H * signal_decay
  # simple 4-neighbor diffusion
  Hnew <- H
  for(i in 2:(G-1)){
    for(j in 2:(G-1)){
      lap <- (H[i-1,j] + H[i+1,j] + H[i,j-1] + H[i,j+1] - 4*H[i,j])
      Hnew[i,j] <- H[i,j] + signal_diffuse * lap
    }
  }
  Hnew
}

## 7) Movement rule (directional bias for E, random for R)
move_agent <- function(x,y, state, H){
  if(state == "D") return(c(x,y))
  p_move <- if(state == "E") p_move_E else p_move_R
  if(runif(1) > p_move) return(c(x,y))
  
  nb <- nbr4(x,y)
  nb[,1] <- clip(nb[,1]); nb[,2] <- clip(nb[,2])
  
  # choose move:
  if(state == "E"){
    # move up gradient of signal (chemotaxis-like)
    vals <- mapply(function(xx,yy) H[xx+1, yy+1], nb[,1], nb[,2])
    nb[which.max(vals), ]
  } else {
    # random walk
    nb[sample(1:4,1), ]
  }
}

## 8) Simulation loop
log_df <- data.frame(step=integer(), nR=integer(), nE=integer(), nD=integer())

for(t in 1:Tmax){
  
  # a) alive agents snapshot
  alive <- agents[agents$state != "D", ]
  
  # b) per-agent contact + alkalinization + signal emission
  for(i in 1:nrow(alive)){
    id <- alive$id[i]
    ax <- alive$x[i]; ay <- alive$y[i]
    
    in_contact <- has_contact(ax, ay, alive[alive$id != id, ])
    edge_alk <- draw_edge_alk(in_contact)
    
    # emit heterogeneity signal if edge alkalinization + contact
    if(edge_alk && in_contact && runif(1) < p_emit_if_edge_alk_contact){
      H[ax+1, ay+1] <- H[ax+1, ay+1] + 0.6
    }
  }
  
  # c) update signal field (decay + diffuse)
  H <- update_signal_field(H)
  
  # d) phenotype transitions + death
  for(i in 1:N){
    if(agents$state[i] == "D") next
    
    ax <- agents$x[i]; ay <- agents$y[i]
    localH <- H[ax+1, ay+1]
    
    # death (can be made condition-dependent)
    if(runif(1) < D_prob_base){
      agents$state[i] <- "D"
      next
    }
    
    # R -> E if signal threshold reached
    if(agents$state[i] == "R" && localH >= E_threshold){
      agents$state[i] <- "E"
    }
  }
  
  # e) movement (avoid collisions: simple sequential update)
  for(i in 1:N){
    if(agents$state[i] == "D") next
    newxy <- move_agent(agents$x[i], agents$y[i], agents$state[i], H)
    
    # collision check: keep if target empty (by alive cells)
    occ <- agents$state != "D"
    if(!any(agents$x[occ] == newxy[1] & agents$y[occ] == newxy[2])){
      agents$x[i] <- newxy[1]
      agents$y[i] <- newxy[2]
    }
  }
  
  # f) log
  log_df <- rbind(log_df, data.frame(
    step=t,
    nR=sum(agents$state=="R"),
    nE=sum(agents$state=="E"),
    nD=sum(agents$state=="D")
  ))
}

## 9) Minimal outputs
print(tail(log_df, 5))

plot(log_df$step, log_df$nE, type="l", xlab="step", ylab="# elongated (E)")

# Optional: visualize final signal field
image(t(H[nrow(H):1, ]), axes=FALSE, main="Final heterogeneity signal field")
box()


#####1 cell
## 0) MINIMAL grid-based ABM (single-cell + population) — sequential script
set.seed(1)

## 1) Parameters
G <- 80                 # tissue grid size (GxG)
N <- 60                 # number of cells (agents)
Tmax <- 400             # steps

p_move_E <- 0.55        # move prob if elongated
p_move_R <- 0.10        # move prob if rounded

# "edge" inside a cell (how many tiles from cell center count as "edge")
cell_radius <- 3
edge_band <- 1

# alkalinization baseline + boosts
p_alk_base <- 0.01
boost_edge <- 0.06
boost_contact <- 0.08

# heterogeneity signal creation
p_emit_if_edge_alk_contact <- 0.35
signal_decay <- 0.92
signal_diffuse <- 0.15

# phenotype transitions
E_threshold <- 1.2      # if local signal exceeds -> R->E
D_prob_base <- 0.001

## 2) Helpers: bounds + neighborhood
clip <- function(x, lo=0, hi=G-1) pmax(lo, pmin(hi, x))

nbr4 <- function(x,y){
  rbind(c(x-1,y), c(x+1,y), c(x,y-1), c(x,y+1))
}

## 3) Initialize population (agent positions, states)
agents <- data.frame(
  id = 1:N,
  x = sample(0:(G-1), N, TRUE),
  y = sample(0:(G-1), N, TRUE),
  state = rep("R", N),          # R, E, D
  stringsAsFactors = FALSE
)

# tissue-level signal field H(x,y)
H <- matrix(0, nrow=G, ncol=G)

## 4) Contact map (cell-cell contact if any other cell in 4-neighborhood)
has_contact <- function(ax, ay, agents_alive){
  nb <- nbr4(ax, ay)
  for(k in 1:nrow(nb)){
    if(any(agents_alive$x == nb[k,1] & agents_alive$y == nb[k,2])) return(TRUE)
  }
  FALSE
}

## 5) Single-cell internal "mini-grid" (lysosome edge alkalinization proxy)
# We don't simulate each lysosome; we simulate a probability summary:
# edge_alk = TRUE/FALSE each time step per cell, driven by edge bias + contact bias.
draw_edge_alk <- function(in_contact){
  p <- p_alk_base + boost_edge + if(in_contact) boost_contact else 0
  runif(1) < pmin(p, 0.95)
}

## 6) Signal diffusion + decay on tissue grid
update_signal_field <- function(H){
  H <- H * signal_decay
  # simple 4-neighbor diffusion
  Hnew <- H
  for(i in 2:(G-1)){
    for(j in 2:(G-1)){
      lap <- (H[i-1,j] + H[i+1,j] + H[i,j-1] + H[i,j+1] - 4*H[i,j])
      Hnew[i,j] <- H[i,j] + signal_diffuse * lap
    }
  }
  Hnew
}

## 7) Movement rule (directional bias for E, random for R)
move_agent <- function(x,y, state, H){
  if(state == "D") return(c(x,y))
  p_move <- if(state == "E") p_move_E else p_move_R
  if(runif(1) > p_move) return(c(x,y))
  
  nb <- nbr4(x,y)
  nb[,1] <- clip(nb[,1]); nb[,2] <- clip(nb[,2])
  
  # choose move:
  if(state == "E"){
    # move up gradient of signal (chemotaxis-like)
    vals <- mapply(function(xx,yy) H[xx+1, yy+1], nb[,1], nb[,2])
    nb[which.max(vals), ]
  } else {
    # random walk
    nb[sample(1:4,1), ]
  }
}

## 8) Simulation loop
log_df <- data.frame(step=integer(), nR=integer(), nE=integer(), nD=integer())

for(t in 1:Tmax){
  
  # a) alive agents snapshot
  alive <- agents[agents$state != "D", ]
  
  # b) per-agent contact + alkalinization + signal emission
  for(i in 1:nrow(alive)){
    id <- alive$id[i]
    ax <- alive$x[i]; ay <- alive$y[i]
    
    in_contact <- has_contact(ax, ay, alive[alive$id != id, ])
    edge_alk <- draw_edge_alk(in_contact)
    
    # emit heterogeneity signal if edge alkalinization + contact
    if(edge_alk && in_contact && runif(1) < p_emit_if_edge_alk_contact){
      H[ax+1, ay+1] <- H[ax+1, ay+1] + 0.6
    }
  }
  
  # c) update signal field (decay + diffuse)
  H <- update_signal_field(H)
  
  # d) phenotype transitions + death
  for(i in 1:N){
    if(agents$state[i] == "D") next
    
    ax <- agents$x[i]; ay <- agents$y[i]
    localH <- H[ax+1, ay+1]
    
    # death (can be made condition-dependent)
    if(runif(1) < D_prob_base){
      agents$state[i] <- "D"
      next
    }
    
    # R -> E if signal threshold reached
    if(agents$state[i] == "R" && localH >= E_threshold){
      agents$state[i] <- "E"
    }
  }
  
  # e) movement (avoid collisions: simple sequential update)
  for(i in 1:N){
    if(agents$state[i] == "D") next
    newxy <- move_agent(agents$x[i], agents$y[i], agents$state[i], H)
    
    # collision check: keep if target empty (by alive cells)
    occ <- agents$state != "D"
    if(!any(agents$x[occ] == newxy[1] & agents$y[occ] == newxy[2])){
      agents$x[i] <- newxy[1]
      agents$y[i] <- newxy[2]
    }
  }
  
  # f) log
  log_df <- rbind(log_df, data.frame(
    step=t,
    nR=sum(agents$state=="R"),
    nE=sum(agents$state=="E"),
    nD=sum(agents$state=="D")
  ))
}

## 9) Minimal outputs
print(tail(log_df, 5))

plot(log_df$step, log_df$nE, type="l", xlab="step", ylab="# elongated (E)")

# Optional: visualize final signal field
image(t(H[nrow(H):1, ]), axes=FALSE, main="Final heterogeneity signal field")
box()

### baska
## Grid-based ABM — full script
## Goal: each cell checks its 8-neighborhood; if ~4 neighbors are occupied (contact)
## and ~4 are empty (free space), it emits a "release/heterogeneity" signal preferentially
## toward the empty side.

set.seed(3)

### 1) Parameters
G <- 90          # grid size (GxG)
N <- 160         # number of cells
Tmax <- 500      # simulation steps

# movement
p_move_R <- 0.10
p_move_E <- 0.55

# signal field dynamics
signal_decay   <- 0.96
signal_diffuse <- 0.10

# release (heterogeneity) signal
release_amount <- 0.8
p_release_base <- 0.05     # baseline when condition is met
p_release_slope <- 0.07    # increases with empty-neighbor count beyond 4

# phenotype transitions
E_threshold <- 1.2
D_prob_base <- 0.001

### 2) Helpers
clip <- function(x, lo=0, hi=G-1) pmax(lo, pmin(hi, x))

nbr8 <- function(x,y){
  rbind(
    c(x-1,y), c(x+1,y), c(x,y-1), c(x,y+1),
    c(x-1,y-1), c(x-1,y+1), c(x+1,y-1), c(x+1,y+1)
  )
}

### 3) Initialize agents
agents <- data.frame(
  id = 1:N,
  x = sample(0:(G-1), N, TRUE),
  y = sample(0:(G-1), N, TRUE),
  state = rep("R", N),   # R, E, D
  stringsAsFactors = FALSE
)

# signal field H(x,y)
H <- matrix(0, nrow=G, ncol=G)

### 4) Occupancy utility
is_occupied <- function(xx, yy, others){
  any(others$x == xx & others$y == yy)
}

### 5) Count empty/occupied neighbors and find an "empty direction"
neighbor_stats <- function(ax, ay, others){
  nb <- nbr8(ax, ay)
  nb[,1] <- clip(nb[,1]); nb[,2] <- clip(nb[,2])
  
  occ <- apply(nb, 1, function(v) is_occupied(v[1], v[2], others))
  empty <- !occ
  
  list(
    nb = nb,
    empty_idx = which(empty),
    occ_idx   = which(occ),
    n_empty = sum(empty),
    n_occ   = sum(occ)
  )
}

# choose one empty neighbor location to deposit release signal (toward free space)
pick_empty_target <- function(stats){
  if(length(stats$empty_idx) == 0) return(NULL)
  stats$nb[sample(stats$empty_idx, 1), , drop=FALSE]
}

### 6) Signal update (decay + diffusion)
update_signal_field <- function(H){
  H <- H * signal_decay
  Hnew <- H
  for(i in 2:(G-1)){
    for(j in 2:(G-1)){
      lap <- (H[i-1,j] + H[i+1,j] + H[i,j-1] + H[i,j+1] - 4*H[i,j])
      Hnew[i,j] <- H[i,j] + signal_diffuse * lap
    }
  }
  Hnew
}

### 7) Movement rule
nbr4 <- function(x,y) rbind(c(x-1,y), c(x+1,y), c(x,y-1), c(x,y+1))

move_agent <- function(x,y,state,H){
  if(state=="D") return(c(x,y))
  p_move <- if(state=="E") p_move_E else p_move_R
  if(runif(1) > p_move) return(c(x,y))
  
  nb <- nbr4(x,y)
  nb[,1] <- clip(nb[,1]); nb[,2] <- clip(nb[,2])
  
  if(state=="E"){
    vals <- mapply(function(xx,yy) H[xx+1, yy+1], nb[,1], nb[,2])
    nb[which.max(vals), ]
  } else {
    nb[sample(1:4, 1), ]
  }
}

### 8) Main loop
log_df <- data.frame(step=integer(), nR=integer(), nE=integer(), nD=integer(),
                     releases=integer())

for(t in 1:Tmax){
  
  alive <- agents[agents$state!="D", ]
  
  ## A) RELEASE events driven by 4-contact / 4-free-space idea
  release_count <- 0
  
  for(i in 1:nrow(alive)){
    id <- alive$id[i]
    ax <- alive$x[i]; ay <- alive$y[i]
    others <- alive[alive$id != id, ]
    
    st <- neighbor_stats(ax, ay, others)
    
    # "4 with contact, 4 empty" approx: require >=4 occupied AND >=4 empty
    if(st$n_occ >= 4 && st$n_empty >= 4){
      
      # probability increases with how "open" the neighborhood is beyond 4 empties
      p_rel <- p_release_base + p_release_slope * (st$n_empty - 4)
      p_rel <- pmin(p_rel, 0.95)
      
      if(runif(1) < p_rel){
        target <- pick_empty_target(st)
        if(!is.null(target)){
          tx <- target[1,1]; ty <- target[1,2]
          H[tx+1, ty+1] <- H[tx+1, ty+1] + release_amount
          release_count <- release_count + 1
        }
      }
    }
  }
  
  ## B) Update signal field
  H <- update_signal_field(H)
  
  ## C) Phenotype transitions + death
  for(i in 1:N){
    if(agents$state[i]=="D") next
    ax <- agents$x[i]; ay <- agents$y[i]
    localH <- H[ax+1, ay+1]
    
    if(runif(1) < D_prob_base){
      agents$state[i] <- "D"
      next
    }
    if(agents$state[i]=="R" && localH >= E_threshold){
      agents$state[i] <- "E"
    }
  }
  
  ## D) Movement with collision avoidance
  for(i in 1:N){
    if(agents$state[i]=="D") next
    newxy <- move_agent(agents$x[i], agents$y[i], agents$state[i], H)
    
    occ <- agents$state!="D"
    if(!any(agents$x[occ]==newxy[1] & agents$y[occ]==newxy[2])){
      agents$x[i] <- newxy[1]
      agents$y[i] <- newxy[2]
    }
  }
  
  ## E) Log
  log_df <- rbind(log_df, data.frame(
    step=t,
    nR=sum(agents$state=="R"),
    nE=sum(agents$state=="E"),
    nD=sum(agents$state=="D"),
    releases=release_count
  ))
}

### 9) Outputs
print(tail(log_df, 6))

plot(log_df$step, log_df$releases, type="l", xlab="step", ylab="# release events")

plot(log_df$step, log_df$nE, type="l", xlab="step", ylab="# elongated (E)")

image(t(H[nrow(H):1, ]), axes=FALSE, main="Final heterogeneity/release signal field")
box()



###
### MULTI-CELL ABM — MORE RED HOTSPOTS VERSION
set.seed(2)

## 1) Parameters
G <- 90
N <- 220
Tmax <- 500

p_move_E <- 0.55
p_move_R <- 0.10

boundary_size <- 12

# ↑ daha fazla olay
p_alk_base <- 0.02
boost_edge <- 0.10
boost_contact <- 0.15

# ↑ emission sıklığı
p_emit_if_edge_alk_contact <- 0.70

# ↓ sönüm, ↓ yayılma → lokal kırmızı
signal_decay <- 0.96
signal_diffuse <- 0.07

# elongation kolaylaşsın
E_threshold <- 0.8
D_prob_base <- 0.001

## 2) Helpers
clip <- function(x, lo=0, hi=G-1) pmax(lo, pmin(hi, x))
nbr4 <- function(x,y) rbind(c(x-1,y), c(x+1,y), c(x,y-1), c(x,y+1))

is_peripheral <- function(x, y, G, b) {
  (x < b) || (x >= G - b) || (y < b) || (y >= G - b)
}

## 3) Init agents
agents <- data.frame(
  id = 1:N,
  x = sample(0:(G-1), N, TRUE),
  y = sample(0:(G-1), N, TRUE),
  state = rep("R", N),
  stringsAsFactors = FALSE
)

H <- matrix(0, nrow=G, ncol=G)

## 4) Contact check
has_contact <- function(ax, ay, others){
  nb <- nbr4(ax, ay)
  for(k in 1:nrow(nb)){
    if(any(others$x == nb[k,1] & others$y == nb[k,2])) return(TRUE)
  }
  FALSE
}

## 5) Edge-alk draw
draw_edge_alk <- function(in_contact){
  p <- p_alk_base + boost_edge + if(in_contact) boost_contact else 0
  runif(1) < pmin(p, 0.98)
}

## 6) Signal update
update_signal_field <- function(H){
  H <- H * signal_decay
  Hnew <- H
  for(i in 2:(G-1)){
    for(j in 2:(G-1)){
      lap <- (H[i-1,j] + H[i+1,j] + H[i,j-1] + H[i,j+1] - 4*H[i,j])
      Hnew[i,j] <- H[i,j] + signal_diffuse * lap
    }
  }
  Hnew
}

## 7) Movement
move_agent <- function(x,y,state,H){
  if(state=="D") return(c(x,y))
  p_move <- if(state=="E") p_move_E else p_move_R
  if(runif(1) > p_move) return(c(x,y))
  
  nb <- nbr4(x,y)
  nb[,1] <- clip(nb[,1]); nb[,2] <- clip(nb[,2])
  
  if(state=="E"){
    vals <- mapply(function(xx,yy) H[xx+1,yy+1], nb[,1], nb[,2])
    nb[which.max(vals), ]
  } else {
    nb[sample(1:4,1), ]
  }
}

## 8) Simulation
log_df <- data.frame(step=integer(), nR=integer(), nE=integer(), nD=integer())

for(t in 1:Tmax){
  
  alive <- agents[agents$state!="D", ]
  
  ## EMISSION
  for(i in 1:nrow(alive)){
    ax <- alive$x[i]; ay <- alive$y[i]
    others <- alive[-i, ]
    
    in_contact <- has_contact(ax, ay, others)
    edge_alk <- draw_edge_alk(in_contact)
    periph <- is_peripheral(ax, ay, G, boundary_size)
    
    if(edge_alk && in_contact && runif(1) < p_emit_if_edge_alk_contact){
      H[ax+1, ay+1] <- H[ax+1, ay+1] + 0.9   # 🔥 daha güçlü sinyal
      if(periph){
        H[ax+1, ay+1] <- H[ax+1, ay+1] + 0.3
      }
    }
  }
  
  H <- update_signal_field(H)
  
  ## STATE TRANSITIONS
  for(i in 1:N){
    if(agents$state[i]=="D") next
    ax <- agents$x[i]; ay <- agents$y[i]
    localH <- H[ax+1, ay+1]
    
    if(runif(1) < D_prob_base){
      agents$state[i] <- "D"; next
    }
    if(agents$state[i]=="R" && localH >= E_threshold){
      agents$state[i] <- "E"
    }
  }
  
  ## MOVE
  for(i in 1:N){
    if(agents$state[i]=="D") next
    newxy <- move_agent(agents$x[i], agents$y[i], agents$state[i], H)
    occ <- agents$state!="D"
    if(!any(agents$x[occ]==newxy[1] & agents$y[occ]==newxy[2])){
      agents$x[i] <- newxy[1]
      agents$y[i] <- newxy[2]
    }
  }
  
  log_df <- rbind(log_df, data.frame(
    step=t,
    nR=sum(agents$state=="R"),
    nE=sum(agents$state=="E"),
    nD=sum(agents$state=="D")
  ))
}

## 9) Plots
plot(log_df$step, log_df$nE, type="l", xlab="step", ylab="# elongated (E)")

image(
  t(H[nrow(H):1, ]),
  axes=FALSE,
  col = colorRampPalette(c("lightyellow","orange","red","darkred"))(50),
  main="Final heterogeneity signal field — enhanced hotspots"
)
box()


