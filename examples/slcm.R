# Standard LCA
slcm(L[3] ~ y1 + y2 + y3)
# Latent transition analysis
slcm(L1[3] ~ y11 + y21 + y31, L2[3] ~ y12 + y22 + y32,
     L1 ~ L2)
# LTA with measurement assumption
slcm(L1[3] ~ y11 + y21 + y31, L2[3] ~ y12 + y22 + y32,
     L1 ~ L2, constraints = c("L1", "L2"))
# Joint latent class analysis
slcm(LX[3] ~ x1 + x2 + x3, LY[3] ~ y1 + y2 + y3,
     LZ[3] ~ z1 + z2 + z3, JL ~ LX + LY + LZ)
