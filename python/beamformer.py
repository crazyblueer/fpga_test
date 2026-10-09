"""
MỤC TIÊU: xuất ra MỘT bộ trọng số hybrid (F_RF, w_BB) thực thi được trên
phần cứng 8-RF, sao cho búp sóng HƯỚNG tới user và đặt NULL vào jammer.

Luồng:
  1. Tái tạo covariance 64 chiều R64 (quét chùm).
  2. Pilot -> góc user;  eigenvector trội -> hướng jammer (để vẽ/kiểm chứng).
  3. Thiết kế trọng số 64 chiều w64 = LCMV:  giữ user, triệt jammer.
  4. GẬP w64 về hybrid (F_RF phase-only khối chéo + w_BB 8 chiều).
  5. Kiểm chứng: vẽ beampattern -> null có nằm đúng chỗ jammer không.

Output cuối: F_RF, w_BB  (và w_hybrid = F_RF @ w_BB để dùng/đánh giá).
"""

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(3)

fc = 28e9; c = 3e8; lam = c / fc; d = lam / 2
N = 64; N_RF = 8; D = N // N_RF
snapshots = 3000
noise_power = 1.0
SNR_dB = 0; JSR_dB = 20
signal_power = 10 ** (SNR_dB / 10)
jammer_power = 10 ** (JSR_dB / 10) * signal_power


def sv(N, ang):
    th = np.deg2rad(ang); n = np.arange(N)
    return (np.exp(1j * 2 * np.pi * d / lam * n * np.sin(th)) / np.sqrt(N)).reshape(-1, 1)


# Kịch bản: jammer gần user để thấy rõ null 
theta_s = 10
theta_jammers = [4, 16]
K = len(theta_jammers)
a_s = sv(N, theta_s)
A_j = [sv(N, t) for t in theta_jammers]

prng = np.random.RandomState(999)
s_pilot = (2 * (prng.randint(0, 2, snapshots) - 0.5) +
           1j * 2 * (prng.randint(0, 2, snapshots) - 0.5)) / np.sqrt(2)
s_sig = np.sqrt(signal_power) * s_pilot.reshape(1, -1)
J = [np.sqrt(jammer_power) *
     ((np.random.randn(snapshots) + 1j * np.random.randn(snapshots)) / np.sqrt(2)).reshape(1, -1)
     for _ in range(K)]
noise = np.sqrt(noise_power / 2) * (np.random.randn(N, snapshots) + 1j * np.random.randn(N, snapshots))
X = a_s @ s_sig
for k in range(K):
    X += A_j[k] @ J[k]
X += noise

# Combiner acquisition -> Y (cho pilot)
beam_angles = np.linspace(-80, 80, N_RF)
F_acq = np.zeros((N, N_RF), dtype=complex)
for r in range(N_RF):
    st, en = r * D, (r + 1) * D
    th = np.deg2rad(beam_angles[r]); nl = np.arange(D)
    F_acq[st:en, r] = np.exp(-1j * 2 * np.pi * d / lam * nl * np.sin(th)) / np.sqrt(D)
Y = F_acq.conj().T @ X


# ----- Tái tạo R64 -----
Q = 2 * D - 1
theta_sweep = np.rad2deg(np.arcsin(-1 + 2 * np.arange(Q) / Q))
E = np.zeros((D * D, 2 * D - 1), dtype=complex)
for i in range(2 * D - 1):
    E[:, i] = np.eye(D, k=i - (D - 1)).flatten()

R64 = np.zeros((N, N), dtype=complex)
for n1 in range(N_RF):
    for n2 in range(N_RF):
        Xn1 = X[n1 * D:(n1 + 1) * D, :]; Xn2 = X[n2 * D:(n2 + 1) * D, :]
        Bm = np.zeros((Q, 2 * D - 1), dtype=complex); pv = np.zeros(Q, dtype=complex)
        for q in range(Q):
            aq = sv(D, theta_sweep[q])
            Bm[q, :] = (np.kron(aq, aq.conj()).T @ E)
            pv[q] = np.mean((aq.conj().T @ Xn1) * (aq.conj().T @ Xn2).conj())
        gamma = np.linalg.pinv(Bm.conj().T @ Bm) @ Bm.conj().T @ pv
        R64[n1 * D:(n1 + 1) * D, n2 * D:(n2 + 1) * D] = (E @ gamma).reshape((D, D), order='F')


# ----- Pilot -> góc user -----
h_pilot = (Y @ s_pilot.conj().reshape(-1, 1)) / snapshots
scan = np.linspace(-90, 90, 2000)
def pilot_score(ang):
    ab = F_acq.conj().T @ sv(N, ang); ab = ab / (np.linalg.norm(ab) + 1e-12)
    return np.abs(ab.conj().T @ h_pilot)[0, 0]
theta_user_hat = scan[np.argmax([pilot_score(a) for a in scan])]
print(f"User thật {theta_s}° -> pilot {theta_user_hat:.1f}° | jammer thật {theta_jammers}")


# ----- Trọng số 64 chiều: LCMV (giữ user, triệt jammer từ R64) -----
# U_d nới rộng quanh góc user để chịu sai số.
U_d = np.hstack([sv(N, theta_user_hat + off) for off in (-2, 0, 2)])
U_d, _ = np.linalg.qr(U_d)
R_inv = np.linalg.inv(R64 + 1e-2 * np.eye(N))
f = np.zeros((U_d.shape[1], 1)); f[1, 0] = 1.0
w64 = R_inv @ U_d @ np.linalg.inv(U_d.conj().T @ R_inv @ U_d) @ f


# ===== GẬP w64
# SPLIT the ideal 64-element weight w64 into a HYBRID pair (F_RF, w_BB)
# that the 8-RF hardware can actually run.
#
# Why we need this:
#   w64 is a perfect 64-element beamformer, but the hardware only has
#   8 RF chains. So we must express it as two stages:
#     F_RF  = analog stage (8 beams, PHASE-ONLY, fixed amplitude)
#     w_BB  = digital stage (8 free complex numbers)
#   Final weight applied by hardware = F_RF @ w_BB.
#
# Key idea (this is what makes it work well):
#   We do NOT try to copy w64. Copying fails because F_RF can only change
#   PHASE, not amplitude, so an exact copy is impossible and the null gets
#   shallow. Instead we RE-SOLVE the real goal directly:
#       "keep the user response = 1, and minimize leftover power
#        (jammer + noise)"  --> this is exactly the MVDR beamformer.
#   We just solve MVDR in the reduced 8-D space that lives AFTER F_RF.
#
# We alternate two easy sub-steps until they settle (15 iterations):
#   (A) fix w_BB  -> best F_RF is just "take the phase" of the weight
#   (B) fix F_RF  -> best w_BB is the MVDR formula (closed form, no search)

def build_FRF_from_phase(w_vec):
    """
    Build F_RF from a 64-element weight by keeping ONLY its phase.
    For each sub-array (each RF chain), take the phase of the matching
    8-element chunk of w_vec; force amplitude to the fixed value 1/sqrt(D)
    so it obeys the phase-shifter (constant-modulus) constraint.
    """
    
    F = np.zeros((N, N_RF), dtype=complex)
    for r in range(N_RF):
        st, en = r * D, (r + 1) * D
        F[st:en, r] = np.exp(1j * np.angle(w_vec[st:en, 0])) / np.sqrt(D)
    return F

F_RF = build_FRF_from_phase(w64)

for _ in range(15):
    # Kênh hiệu dụng sau F_RF: mọi thứ chiếu xuống 8 chiều.
    R_eff = F_RF.conj().T @ R64 @ F_RF + 1e-3 * np.eye(N_RF)   # covariance 8 chiều hiệu dụng
    a_eff = F_RF.conj().T @ sv(N, theta_user_hat)              # user trong 8 chiều hiệu dụng
    # w_BB = MVDR: min công suất, giữ đáp ứng user = 1
    Rinv_eff = np.linalg.inv(R_eff)
    w_BB = Rinv_eff @ a_eff / (a_eff.conj().T @ Rinv_eff @ a_eff)
    # Cập nhật F_RF từ pha của trọng số hybrid hiện tại (tinh chỉnh nhẹ)
    w_now = F_RF @ w_BB
    F_RF = build_FRF_from_phase(w_now)

w_hybrid = F_RF @ w_BB          # <<< BỘ TRỌNG SỐ CUỐI CÙNG (thực thi được)


# ----- Đánh giá -----
def sinr_full(w):
    w = w.reshape(-1, 1)
    s = signal_power * np.abs(w.conj().T @ a_s)[0, 0] ** 2
    it = sum(jammer_power * np.abs(w.conj().T @ A_j[k])[0, 0] ** 2 for k in range(K))
    nz = noise_power * np.real(w.conj().T @ w)[0, 0]
    return 10 * np.log10(s / (it + nz))

w_naive = sv(N, theta_user_hat)
print(f"\nSINR beamformer thường (không triệt): {sinr_full(w_naive):6.2f} dB")
print(f"SINR HYBRID (F_RF, w_BB) cuối cùng:   {sinr_full(w_hybrid):6.2f} dB")
print(f"Cải thiện:                            {sinr_full(w_hybrid) - sinr_full(w_naive):6.2f} dB")


# ----- Kiểm chứng bằng beampattern: null có đúng chỗ jammer không -----
angles = np.linspace(-90, 90, 3000)
pattern = np.array([np.abs(w_hybrid.conj().T @ sv(N, a))[0, 0] for a in angles])
pattern_db = 20 * np.log10(pattern / pattern.max() + 1e-12)

plt.figure(figsize=(11, 5))
plt.plot(angles, pattern_db, 'b', lw=1.3, label='Búp sóng hybrid (F_RF, w_BB)')
plt.axvline(theta_s, color='green', ls='--', lw=2, label=f'User ({theta_s}°)')
for jt in theta_jammers:
    plt.axvline(jt, color='red', ls=':', lw=2, label='Jammer' if jt == theta_jammers[0] else None)
plt.xlabel('Góc (độ)'); plt.ylabel('Đáp ứng búp sóng (dB)')
plt.title('Búp sóng: hướng tới user, đặt null vào jammer')
plt.ylim([-60, 3]); plt.xlim([-90, 90]); plt.grid(True, ls='--', alpha=0.5)
plt.legend(loc='lower right'); plt.tight_layout()
#plt.savefig('path to\\beampattern.png', dpi=110)

# In đáp ứng tại các hướng quan tâm để kiểm chứng bằng số
resp_user = 20 * np.log10(np.abs(w_hybrid.conj().T @ sv(N, theta_s))[0, 0] / pattern.max())
print(f"\nĐáp ứng tại user  ({theta_s}°): {resp_user:6.1f} dB  (càng gần 0 càng tốt)")
for jt in theta_jammers:
    r = 20 * np.log10(np.abs(w_hybrid.conj().T @ sv(N, jt))[0, 0] / pattern.max())
    print(f"Đáp ứng tại jammer ({jt}°): {r:6.1f} dB  (càng thấp = null càng sâu)")

# Xuất bộ trọng số ra file để dùng/nạp lên phần cứng
#np.savez('path to\\weights.npz', F_RF=F_RF, w_BB=w_BB, w_hybrid=w_hybrid)
print("\nĐã lưu bộ trọng số: weights.npz  (F_RF, w_BB, w_hybrid)")
#print("Đã lưu hình búp sóng: path to\\beampattern.png")

####

# ===== FIXED-POINT MODEL + TEST VECTORS =====
def q15(v): return np.clip(np.round(v * 32768), -32768, 32767).astype(np.int64)
PH = 6                      # phase bits
M  = 200                    # number of test snapshots

# --- coefficients (already conjugated: hardware computes y = sum(c*x)) ---
cidx = np.zeros(N, dtype=int)
for r in range(N_RF):
    for k in range(D):
        n = r * D + k
        cidx[n] = int(np.round(np.angle(F_RF[n, r]) / (2*np.pi / 2**PH))) % 2**PH
c  = np.exp(-1j * 2*np.pi * cidx / 2**PH)
b  = np.conj(w_BB.flatten()) / np.sqrt(D)
b  = b * 0.99 / np.max(np.abs(b))
cr, ci = q15(c.real), q15(c.imag)
br, bi = q15(b.real), q15(b.imag)

# --- check 1: do the quantized weights keep the nulls? ---
v = np.array([(br[n//D] + 1j*bi[n//D]) * (cr[n] + 1j*ci[n]) for n in range(N)]) / 32768**2
def resp(a): return np.abs(v @ sv(N, a)[:, 0])
ref = resp(theta_s)
print("\n--- quantized weights (6-bit phase, Q1.15) ---")
print(f"user   {theta_s}°: 0.0 dB")
for jt in theta_jammers:
    print(f"jammer {jt}°: {20*np.log10(resp(jt)/ref):6.1f} dB")

# --- check 2: bit-exact golden outputs ---
xr, xi = q15(X.real[:, :M] / 16), q15(X.imag[:, :M] / 16)
yg = []
for t in range(M):
    zr = np.zeros(N_RF, dtype=np.int64); zi = np.zeros(N_RF, dtype=np.int64)
    for r in range(N_RF):
        s = slice(r*D, (r+1)*D)
        zr[r] = np.sum(cr[s]*xr[s, t] - ci[s]*xi[s, t]) >> 15
        zi[r] = np.sum(cr[s]*xi[s, t] + ci[s]*xr[s, t]) >> 15
    yr = int(np.sum(br*zr - bi*zi)) >> 15
    yi = int(np.sum(br*zi + bi*zr)) >> 15
    yg.append((yr, yi))

h16 = lambda v: f"{int(v) & 0xFFFF:04x}"
with open('x_vec.hex', 'w') as f:
    for t in range(M):
        for n in range(N): f.write(h16(xr[n, t]) + h16(xi[n, t]) + "\n")
with open('coef.hex', 'w') as f:
    for n in range(N): f.write(h16(cr[n]) + h16(ci[n]) + "\n")
with open('bb.hex', 'w') as f:
    for r in range(N_RF): f.write(h16(br[r]) + h16(bi[r]) + "\n")
with open('y_gold.hex', 'w') as f:
    for yr, yi in yg: f.write(f"{yr & 0xFFFFFF:06x}{yi & 0xFFFFFF:06x}\n")
print("wrote x_vec.hex coef.hex bb.hex y_gold.hex")