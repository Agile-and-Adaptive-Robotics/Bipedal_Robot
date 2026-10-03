L = open(r'D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo\minimizeFlxPin.m', encoding='utf8', errors='replace').read().splitlines()
for a, b in [(415, 426), (434, 441), (458, 466), (475, 480), (538, 546)]:
    print(f'--- minimizeFlxPin {a}-{b} ---')
    print("\n".join(f'{i+1}: {L[i]}' for i in range(a - 1, b)))
M = open(r'D:\Github\Bipedal_Robot\Code\Matlab\Robot_Data\MonoPam_mult.m', encoding='utf8', errors='replace').read().splitlines()
print('--- MonoPam_mult 596-612 ---')
print("\n".join(f'{i+1}: {M[i]}' for i in range(595, 612)))
hits = [(i + 1, l.strip()) for i, l in enumerate(M) if 'fzero' in l.lower()]
print('fzero in MonoPam_mult:', hits[:8])
