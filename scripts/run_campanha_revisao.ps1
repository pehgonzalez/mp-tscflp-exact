#Requires -Version 5.1
<#
.SYNOPSIS
  Campanha unica da revisao EJOR, build fresco, com retomada segura.

.DESCRIPTION
  Cobre TODOS os experimentos pendentes das duas rodadas de parecer:

    bloco fatorial-50-5       R8 + item 2 da rodada 2. Cinco bracos no mesmo
                              build sobre as 25 instancias 50-5 a 1800 s:
                              BBC, L-BBC completo, L-BBC sem warm start,
                              L-BBC sem master start (reproducao do defeito),
                              L-BBC sem a linha do bound. 125 runs.
    bloco 50-10-1800          R10. O conjunto 50-10 no orcamento cheio,
                              3 engines. 75 runs. Tambem elimina os runs de
                              decomposicao de 10/07, anteriores ao build warm.
    bloco 100-sem-corepoint   R9. Os dois conjuntos de 100 plantas com as
                              duas decomposicoes SEM core-point. 100 runs.
    bloco baseline-fl         R13a. MIP compacto com as linhas (F_l). 100 runs.
    bloco baseline-mipguided  R13b. MIP compacto com v_LD, start reparado e
                              cutoff. 100 runs.
    bloco variancia           R12. 12 instancias abertas (PSC1-C2/C4/C5 de
                              cada conjunto), engine de melhor gap, 3
                              repeticoes cada (exact-var1..3). 36 runs.
    bloco protocolo-mauri     R16. As 12 fechadas no protocolo de Mauri et
                              al., 3600 s, tolerancias default do solver,
                              sem proof mode. 12 runs.
    bloco estendido-confirma  R18. Os 2 runs estendidos de 19/07 refeitos no
                              build atual, BBC 14400 s. 2 runs.
    bloco combo-fl-guided     Extensao pos-analise. MIP compacto com (F_l) e
                              v_LD juntos, as 100 instancias a 1800 s.
    bloco proveniencia-cp     Extensao pos-analise. Fecha as 12 sob o build
                              lancado, com core-point ligado como nos
                              fechamentos originais. 3 runs (1800 s + 2x14400 s).

  Blocos da rodada 3 do parecer (2026-08-17), todos depois dos anteriores:

    bloco cert-dump           R2. As 12 fechadas de novo, cada uma com a engine
                              que a fechou no lote canonico, com
                              MPTSCFL_CERT_DUMP apontando para logs\certdump.
                              Rotulo exact-cert3. 12 runs, 8 a 1800 s e as 4
                              que so fecham no orcamento estendido a 14400 s.
                              20.0 h. O env var e setado so neste bloco e
                              limpo depois de cada run.
    bloco mip-fresco          R7. O MIP compacto (m0, exact) nas 100
                              instancias a 1800 s, no build fresco. 100 runs,
                              50.0 h. O conjunto 50-10 ja rodou pos-epoca com
                              esta chave exata e e pulado sozinho, entao na
                              pratica sao 75 runs e 37.5 h.
    bloco fatorial-100        R8. O fatorial de 5 bracos do bloco
                              fatorial-50-5 repetido no conjunto 100-10.
                              125 runs, 62.5 h. Os bracos m1 exact e m2 exact
                              coincidem com o 100-sem-corepoint (mesma chave,
                              mesma configuracao) e sao pulados se aquele
                              bloco ja rodou, o que deixa 75 runs e 37.5 h.
    bloco oneeval             R3. m2 exact-oneeval, uma unica avaliacao
                              Lagrangiana nos duais de LP. 100 runs, 50.0 h.
    bloco repaircb            R5. m0 exact-repaircb, MIP compacto com o
                              callback de reparo por fluxo. 100 runs, 50.0 h.
    bloco thetaint            R12. m2 exact-thetaint, theta_l inteiro no
                              master. 100 runs, 50.0 h. EXIGE o patch
                              patch_thetaint.md aplicado: sem ele o binario
                              sai com codigo 3 e nao grava linha.
    bloco flgcb               m0 exact-flgcb, o combo (F_l)+v_LD com o
                              callback de reparo. 100 runs, 50.0 h.
    bloco sementes            R15. As 12 instancias do bloco variancia, cada
                              uma na engine da tabela VarPlan, modo
                              exact-var1, seeds 1 a 5. 60 runs, 30.0 h.

  1350 runs no plano completo. Os 653 anteriores somam 288 h mais 58.5 h da
  extensao; os 697 da rodada 3 somam 362.5 h de plano. Contra um results.csv
  que ja tenha o lote canonico pos-epoca, 75 desses runs sao pulados sozinhos
  (25 do mip-fresco, do conjunto 50-10, e 50 do fatorial-100, os bracos m1 e
  m2 exact que o bloco 100-sem-corepoint ja rodou), o que deixa 622 runs e
  ~325 h. Runs ja gravados desde a epoca sao pulados automaticamente.

  SEEDS. Ate a rodada 2 todo run usava o -Seed global e o filtro de
  "ja feito" so olhava essa seed. O bloco sementes quebra isso, entao cada
  item do plano carrega o proprio campo Seed (default, o -Seed global) e a
  chave de "ja feito" passa a incluir a seed da LINHA do results.csv. Um run
  de seed 3 nao mascara mais o mesmo run em seed 0 nem vice-versa.

  RETOMADA. A epoca do build fresco fica em logs\campanha_revisao\epoch.txt.
  Um run e pulado quando logs\results.csv ja tem uma linha com a mesma
  instancia, mode, method e time limit datada da epoca em diante. Pode
  interromper com Ctrl-C ou criando o arquivo logs\campanha_revisao\STOP, e
  relancar retoma do ponto exato.

  PROGRESSO. Cada run imprime bloco, posicao na fila, instancia, hora de
  entrada, time limit e a estimativa de termino da campanha inteira,
  refinada pela razao media entre tempo real e time limit dos runs ja
  concluidos nesta sessao.

.PARAMETER Mode
  Report  imprime o plano e o que ja esta feito, nao executa nada
  Apply   executa
.PARAMETER Blocks
  Lista separada por virgula para rodar so alguns blocos. Padrao, todos.
.EXAMPLE
  powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\run_campanha_revisao.ps1 -Mode Report
  powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\run_campanha_revisao.ps1 -Mode Apply
#>
[CmdletBinding()]
param(
    [ValidateSet('Report','Apply')]
    [string]$Mode = 'Report',
    [string]$Blocks = '',
    [int]$Seed = 0,
    [int]$Threads = 16,
    # So e preciso quando a variavel de ambiente GUROBI_HOME nao existe.
    # Exemplo: -GurobiHome C:\gurobi1302\win64
    [string]$GurobiHome = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version 2.0

# O binario grava logs\results.csv relativo ao diretorio corrente, entao o
# script se muda para a raiz do repositorio, a pasta pai de scripts\.
if ($PSScriptRoot) { Set-Location (Split-Path -Parent $PSScriptRoot) }
$Root = (Get-Location).Path

$Exe      = Join-Path $Root 'build\Release\mptscfl.exe'
$TestExe  = Join-Path $Root 'build\Release\test_core.exe'
$DataDir  = 'data\instances'
$Results  = Join-Path $Root 'logs\results.csv'
$StateDir = Join-Path $Root 'logs\campanha_revisao'
$EpochF   = Join-Path $StateDir 'epoch.txt'
$StopF    = Join-Path $StateDir 'STOP'
$Journal  = Join-Path $StateDir 'journal.csv'

function Say([string]$m) { Write-Host ("[{0}] {1}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $m) }

# ---------------------------------------------------------------- o plano ----

$Sets   = @('50-5','50-10','100-5','100-10')
$Closed = @('PSC1-C1-100-5','PSC1-C1-50-5','PSC2-C1-100-5','PSC2-C1-50-5',
            'PSC3-C1-100-10','PSC3-C1-100-5','PSC3-C1-50-5','PSC4-C1-100-10',
            'PSC4-C1-100-5','PSC4-C1-50-5','PSC5-C1-100-5','PSC5-C1-50-5')
# Instancia da variancia -> engine de melhor gap canonico (calculado do
# results.csv em 2026-08-05; ver AVALIACAO/RESPONSE da rodada 2).
$VarPlan = @(
    @{ I='PSC1-C2-50-5';   M=0 }, @{ I='PSC1-C4-50-5';   M=0 }, @{ I='PSC1-C5-50-5';   M=0 },
    @{ I='PSC1-C2-50-10';  M=0 }, @{ I='PSC1-C4-50-10';  M=2 }, @{ I='PSC1-C5-50-10';  M=2 },
    @{ I='PSC1-C2-100-5';  M=0 }, @{ I='PSC1-C4-100-5';  M=2 }, @{ I='PSC1-C5-100-5';  M=2 },
    @{ I='PSC1-C2-100-10'; M=2 }, @{ I='PSC1-C4-100-10'; M=2 }, @{ I='PSC1-C5-100-10'; M=2 })

# Bloco cert-dump. Para cada uma das 12 fechadas, a engine que a fechou no lote
# canonico (a mais rapida entre as que fecharam) e o orcamento em que fechou.
# As 4 ultimas so fecham no orcamento estendido, entao vao a 14400 s.
# Excecao anotada: PSC3-C1-100-10 nao fechou no lote canonico a 1800 s; o unico
# fechamento a 1800 s dela e o run m1 anterior a epoca, o mesmo que o bloco
# proveniencia-cp esta refazendo. Mantida a 1800 s com m1; se o run do
# cert-dump nao fechar, o certificado dela sai do run de 14400 s que for
# lancado a parte, nao deste bloco.
$CertPlan = @(
    @{ I='PSC1-C1-50-5';   M=1; TL=1800  }, @{ I='PSC1-C1-100-5';  M=2; TL=1800  },
    @{ I='PSC3-C1-100-5';  M=2; TL=1800  }, @{ I='PSC3-C1-100-10'; M=1; TL=1800  },
    @{ I='PSC4-C1-50-5';   M=1; TL=1800  }, @{ I='PSC4-C1-100-5';  M=2; TL=1800  },
    @{ I='PSC5-C1-50-5';   M=1; TL=1800  }, @{ I='PSC5-C1-100-5';  M=1; TL=1800  },
    @{ I='PSC2-C1-50-5';   M=1; TL=14400 }, @{ I='PSC2-C1-100-5';  M=2; TL=14400 },
    @{ I='PSC3-C1-50-5';   M=1; TL=14400 }, @{ I='PSC4-C1-100-10'; M=2; TL=14400 })

function All-Instances([string]$suffix) {
    $out = @()
    foreach ($p in 1..5) { foreach ($c in 1..5) { $out += "PSC$p-C$c-$suffix" } }
    return $out
}

function Build-Plan {
    $plan = New-Object System.Collections.ArrayList
    # $sd e a seed do run; default, a -Seed global. So o bloco sementes passa
    # outra coisa. Todo item do plano carrega o campo, entao Key pode le-lo sem
    # tropecar no Set-StrictMode.
    $add = { param($b,$i,$m,$md,$tl,$pap,$sd=$Seed)
        [void]$plan.Add([pscustomobject]@{ Bloco=$b; Inst=$i; Method=$m; Md=$md; TL=$tl; Pap=$pap; Seed=$sd }) }

    foreach ($i in (All-Instances '50-5')) {
        & $add 'fatorial-50-5' $i 1 'exact'         1800 0
        & $add 'fatorial-50-5' $i 2 'exact'         1800 0
        & $add 'fatorial-50-5' $i 2 'exact-nowarm'  1800 0
        & $add 'fatorial-50-5' $i 2 'exact-nostart' 1800 0
        & $add 'fatorial-50-5' $i 2 'exact-noldb'   1800 0
    }
    foreach ($i in (All-Instances '50-10')) {
        foreach ($m in 0,1,2) { & $add '50-10-1800' $i $m 'exact' 1800 0 }
    }
    foreach ($sfx in '100-5','100-10') {
        foreach ($i in (All-Instances $sfx)) {
            foreach ($m in 1,2) { & $add '100-sem-corepoint' $i $m 'exact' 1800 0 }
        }
    }
    foreach ($sfx in $Sets) {
        foreach ($i in (All-Instances $sfx)) { & $add 'baseline-fl' $i 0 'exact-fl' 1800 0 }
    }
    foreach ($sfx in $Sets) {
        foreach ($i in (All-Instances $sfx)) { & $add 'baseline-mipguided' $i 0 'exact-mipguided' 1800 0 }
    }
    foreach ($v in $VarPlan) {
        foreach ($rep in 1,2,3) { & $add 'variancia' $v.I $v.M ("exact-var{0}" -f $rep) 1800 0 }
    }
    foreach ($i in $Closed) { & $add 'protocolo-mauri' $i 0 'exact-mauri' 3600 0 }
    & $add 'estendido-confirma' 'PSC2-C1-50-5' 1 'exact' 14400 0
    & $add 'estendido-confirma' 'PSC3-C1-50-5' 1 'exact' 14400 0
    # Blocos da extensao pos-analise (2026-08-17).
    # combo-fl-guided: o MIP compacto com (F_l) E v_LD juntos, a celula que a
    # analise da campanha apontou como a provavel melhor engine da tabela.
    foreach ($sfx in $Sets) {
        foreach ($i in (All-Instances $sfx)) { & $add 'combo-fl-guided' $i 0 'exact-flguided' 1800 0 }
    }
    # proveniencia-cp: fecha as 12 sob o build lancado sem asterisco. Os tres
    # runs rodam com core-point LIGADO, como os fechamentos originais, e levam
    # o rotulo exact-pap para nao colidir com os runs sem cortes ja gravados.
    & $add 'proveniencia-cp' 'PSC3-C1-100-10' 1 'exact-pap' 1800 1
    & $add 'proveniencia-cp' 'PSC2-C1-50-5'   1 'exact-pap' 14400 1
    & $add 'proveniencia-cp' 'PSC3-C1-50-5'   1 'exact-pap' 14400 1
    # ------------------------------------------- blocos da rodada 3 (R2..R15) ----
    # cert-dump: refaz as 12 fechadas com o despejo do certificado racional
    # ligado. O rotulo exact-cert3 e a terceira geracao do experimento (exact-cert
    # de 07/2026 ja existe no results.csv), e nao colide com nenhum modo especial
    # do main.cpp: qualquer rotulo iniciado por "exact" roda o modelo exato puro.
    foreach ($c in $CertPlan) { & $add 'cert-dump' $c.I $c.M 'exact-cert3' $c.TL 0 }
    # mip-fresco: o MIP compacto nas 100 no build fresco. A chave
    # (inst, exact, 0, 1800) ja existe pos-epoca para o conjunto 50-10 e so para
    # ele, entao esses 25 runs sao pulados pelo mecanismo normal, que e o
    # comportamento desejado: aqueles ja sao frescos.
    foreach ($sfx in $Sets) {
        foreach ($i in (All-Instances $sfx)) { & $add 'mip-fresco' $i 0 'exact' 1800 0 }
    }
    # fatorial-100: os mesmos 5 bracos do fatorial-50-5, agora em 100-10.
    foreach ($i in (All-Instances '100-10')) {
        & $add 'fatorial-100' $i 1 'exact'         1800 0
        & $add 'fatorial-100' $i 2 'exact'         1800 0
        & $add 'fatorial-100' $i 2 'exact-nowarm'  1800 0
        & $add 'fatorial-100' $i 2 'exact-nostart' 1800 0
        & $add 'fatorial-100' $i 2 'exact-noldb'   1800 0
    }
    # oneeval / repaircb / thetaint / flgcb: os quatro modos novos do main.cpp,
    # cada um nas 100 instancias a 1800 s.
    foreach ($sfx in $Sets) {
        foreach ($i in (All-Instances $sfx)) { & $add 'oneeval'  $i 2 'exact-oneeval'  1800 0 }
    }
    foreach ($sfx in $Sets) {
        foreach ($i in (All-Instances $sfx)) { & $add 'repaircb' $i 0 'exact-repaircb' 1800 0 }
    }
    foreach ($sfx in $Sets) {
        foreach ($i in (All-Instances $sfx)) { & $add 'thetaint' $i 2 'exact-thetaint' 1800 0 }
    }
    foreach ($sfx in $Sets) {
        foreach ($i in (All-Instances $sfx)) { & $add 'flgcb'    $i 0 'exact-flgcb'    1800 0 }
    }
    # sementes: as 12 da variancia na engine do VarPlan, modo exact-var1, seeds
    # 1..5. A seed 0 fica de fora de proposito: ela ja esta gravada pelo bloco
    # variancia, com esse mesmo rotulo, e entra na analise como sexta repeticao.
    foreach ($v in $VarPlan) {
        foreach ($sd in 1..5) { & $add 'sementes' $v.I $v.M 'exact-var1' 1800 0 $sd }
    }
    return $plan
}

# ------------------------------------------------------- feito ou pendente ----

function Load-Done([string]$epoch) {
    # Chaves instancia|mode|method|tl|seed com linha datada da epoca em diante.
    # A seed entra na CHAVE em vez de filtrar as linhas: o bloco sementes roda o
    # mesmo (instancia, modo, method, tl) em cinco seeds, e com o filtro antigo
    # (so linhas da -Seed global) as seeds 1..5 nunca seriam vistas como feitas,
    # o que faria a campanha repetir os 60 runs a cada relancamento.
    $done = @{}
    if (-not (Test-Path -LiteralPath $Results)) { return $done }
    foreach ($r in (Import-Csv -LiteralPath $Results)) {
        if ($r.datetime -lt $epoch) { continue }
        $tl = [string][int][double]$r.time_limit_s
        $sd = [string][int][double]$r.seed
        $k = '{0}|{1}|{2}|{3}|{4}' -f $r.instance, $r.mode, $r.method, $tl, $sd
        $done[$k] = $true
    }
    return $done
}

function Key($r) { '{0}|{1}|{2}|{3}|{4}' -f $r.Inst, $r.Md, $r.Method, $r.TL, ([int]$r.Seed) }

# ----------------------------------------------------------------- inicio ----

Write-Host ('=' * 78)
Write-Host ("Campanha unica da revisao EJOR, modo {0}" -f $Mode)
Write-Host ('=' * 78)

if (-not (Test-Path -LiteralPath (Join-Path $Root $DataDir))) {
    throw "data\instances nao encontrada em $Root. Rode da raiz do repositorio."
}
$nInst = @(Get-ChildItem -Path (Join-Path $Root $DataDir) -Filter 'PSC*.txt' -File).Count
if ($nInst -ne 100) { throw "data\instances tem $nInst instancias, esperava 100." }

$plan = Build-Plan
$wanted = @()
if ($Blocks.Trim() -ne '') { $wanted = @($Blocks.Split(',') | ForEach-Object { $_.Trim() }) }
if ($wanted.Count -gt 0) { $plan = @($plan | Where-Object { $wanted -contains $_.Bloco }) }

if (-not (Test-Path -LiteralPath $StateDir)) { New-Item -ItemType Directory -Path $StateDir -Force | Out-Null }

# Epoca. Definida uma unica vez, DEPOIS do primeiro build fresco.
$epoch = ''
if (Test-Path -LiteralPath $EpochF) {
    $epoch = (Get-Content -LiteralPath $EpochF | Select-Object -First 1).Trim()
    Say "epoca da campanha ja definida: $epoch (retomada)"
}

if ($Mode -eq 'Apply') {
    if (Test-Path -LiteralPath $StopF) { Remove-Item -LiteralPath $StopF -Force }
    # O repositorio reorganizado nunca teve build\ (a pasta e ignorada pelo
    # git e nao foi migrada), entao a primeira execucao configura o CMake.
    if (-not (Test-Path -LiteralPath (Join-Path $Root 'build\CMakeCache.txt'))) {
        Say 'build\ inexistente, configurando o CMake pela primeira vez neste layout'
        $cfg = @('-S', '.', '-B', 'build')
        if ($GurobiHome -ne '') {
            $cfg += ('-DGUROBI_HOME=' + $GurobiHome)
        } elseif (-not (Test-Path Env:GUROBI_HOME)) {
            Say 'AVISO: variavel GUROBI_HOME nao definida e -GurobiHome nao passado.'
            Say 'Se o configure nao achar o Gurobi, rode de novo com, por exemplo,'
            Say '  -GurobiHome C:\gurobi1302\win64'
        }
        cmake @cfg
        if ($LASTEXITCODE -ne 0) { throw 'CMAKE CONFIGURE FALHOU, nada foi executado.' }
    }
    Say 'build fresco (cmake --build build --config Release)'
    cmake --build build --config Release
    if ($LASTEXITCODE -ne 0) { throw 'BUILD FALHOU, nada foi executado.' }
    if (-not (Test-Path -LiteralPath $Exe)) {
        throw ("binario ausente apos o build: $Exe . Quase sempre isso quer dizer " +
               "que o CMake nao achou o Gurobi e compilou so o test_core. Apague a " +
               "pasta build e rode de novo com -GurobiHome C:\gurobi1302\win64 " +
               "(ajuste a versao para a instalada).")
    }
    if (Test-Path -LiteralPath $TestExe) {
        Say 'teste de unidade (test_core)'
        & $TestExe
        if ($LASTEXITCODE -ne 0) { throw 'test_core FALHOU, campanha abortada.' }
    } else { Say 'test_core.exe ausente, seguindo sem o teste' }
    if ($epoch -eq '') {
        $epoch = Get-Date -Format 'yyyyMMdd-HHmmss'
        Set-Content -LiteralPath $EpochF -Value $epoch -Encoding ASCII
        Say "epoca do build fresco registrada: $epoch"
    }
    if (-not (Test-Path -LiteralPath $Journal)) {
        Set-Content -LiteralPath $Journal -Value 'seq,entrada,bloco,instancia,method,mode,tl,dur_s,exit,acao' -Encoding UTF8
    }
}
if ($epoch -eq '') { $epoch = '99999999-999999' } # Report antes do 1o Apply: nada feito

$done = Load-Done $epoch
$pend = @($plan | Where-Object { -not $done.ContainsKey((Key $_)) })
$feitos = $plan.Count - $pend.Count

$somaTL = 0.0
foreach ($r in $pend) { $somaTL += [double]$r.TL }

Write-Host ''
Say ("plano total {0} runs | concluidos {1} | pendentes {2}" -f $plan.Count, $feitos, $pend.Count)
foreach ($g in ($plan | Group-Object Bloco)) {
    $pg = @($pend | Where-Object { $_.Bloco -eq $g.Name })
    $h = 0.0; foreach ($r in $pg) { $h += [double]$r.TL }
    Write-Host ("  {0,-20} {1,3} runs, pendentes {2,3}, ~{3,6:N1} h" -f $g.Name, $g.Count, $pg.Count, ($h/3600))
}
Say ("orcamento pendente somado: {0:N1} h (~{1:N1} dias)" -f ($somaTL/3600), ($somaTL/86400))

if ($Mode -eq 'Report') {
    Write-Host ''
    Write-Host 'Nada foi executado. Para rodar:' -ForegroundColor Green
    Write-Host '  powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\run_campanha_revisao.ps1 -Mode Apply' -ForegroundColor Green
    return
}

# ------------------------------------------------------------ o laco real ----

$ratio = 1.0      # razao media tempo real / TL, refinada a cada run
$nMed  = 0
$seq   = 0
$falhas = 0
$campT0 = Get-Date

for ($idx = 0; $idx -lt $pend.Count; $idx++) {
    $r = $pend[$idx]

    if (Test-Path -LiteralPath $StopF) {
        Say 'STOP encontrado, parando com seguranca. Relancar retoma daqui.'
        break
    }
    # Alguem pode ter completado este run em outra janela; reconfere barato.
    $done2 = Load-Done $epoch
    if ($done2.ContainsKey((Key $r))) {
        Say ("pulado (ja no results.csv): {0} {1} m{2}" -f $r.Inst, $r.Md, $r.Method)
        continue
    }

    $restante = 0.0
    for ($j = $idx; $j -lt $pend.Count; $j++) { $restante += [double]$pend[$j].TL }
    $etaS = $restante * $ratio
    $eta  = (Get-Date).AddSeconds($etaS)
    $entrada = Get-Date

    Write-Host ''
    Write-Host ('-' * 78)
    Say ("[{0}] run {1}/{2} da fila | {3} | m{4} {5} | TL {6} s | seed {7}" -f `
        $r.Bloco, ($idx+1), $pend.Count, $r.Inst, $r.Method, $r.Md, $r.TL, $r.Seed)
    Say ("entrada {0:HH:mm:ss} | restante ~{1:N1} h | termino estimado da campanha {2:ddd dd/MM HH:mm}" -f `
        $entrada, ($etaS/3600), $eta)

    # Despejo do certificado racional: ligado SO no bloco cert-dump e limpo
    # logo depois do run, para nenhum outro bloco gravar .certdump por acidente.
    $certOn = ($r.Bloco -eq 'cert-dump')
    if (Test-Path Env:MPTSCFL_CERT_DUMP) { Remove-Item Env:MPTSCFL_CERT_DUMP }
    if ($certOn) {
        $certDir = 'logs\certdump'
        if (-not (Test-Path -LiteralPath (Join-Path $Root $certDir))) {
            New-Item -ItemType Directory -Path (Join-Path $Root $certDir) -Force | Out-Null
        }
        $env:MPTSCFL_CERT_DUMP = $certDir
        Say ("MPTSCFL_CERT_DUMP={0} (so neste bloco)" -f $certDir)
    }

    $instFile = Join-Path $DataDir ($r.Inst + '.txt')
    try {
        & $Exe $instFile $r.Method $r.TL $r.Md $r.Seed $Threads $r.Pap
        $rc = $LASTEXITCODE
    } finally {
        if ($certOn -and (Test-Path Env:MPTSCFL_CERT_DUMP)) { Remove-Item Env:MPTSCFL_CERT_DUMP }
    }
    $dur = ((Get-Date) - $entrada).TotalSeconds

    $seq++
    Add-Content -LiteralPath $Journal -Encoding UTF8 -Value `
        ('{0},{1:yyyy-MM-ddTHH:mm:ss},{2},{3},{4},{5},{6},{7:N0},{8},{9}' -f `
        $seq, $entrada, $r.Bloco, $r.Inst, $r.Method, $r.Md, $r.TL, $dur, $rc,
        ('executado seed=' + $r.Seed))
    # A seed vai dentro da coluna acao de proposito: acrescentar uma coluna
    # nova quebraria os journal.csv ja gravados com o cabecalho antigo.

    if ($rc -ne 0) {
        $falhas++
        Say ("AVISO: exit {0} em {1} m{2} {3}; seguindo (o run pode ser refeito relancando)" -f $rc, $r.Inst, $r.Method, $r.Md)
    }
    # Refina a razao real/TL com media movel (runs fechados cedo puxam para baixo).
    $obs = $dur / [double]$r.TL
    if ($obs -gt 0.01) {
        $nMed++
        $ratio = (($ratio * ($nMed - 1)) + [Math]::Min($obs, 1.2)) / $nMed
        if ($ratio -lt 0.25) { $ratio = 0.25 }
    }
}

Write-Host ''
Write-Host ('=' * 78)
$dTot = ((Get-Date) - $campT0)
Say ("sessao encerrada | runs executados {0} | falhas {1} | duracao {2:N1} h" -f $seq, $falhas, $dTot.TotalHours)
$done = Load-Done $epoch
$pendF = @($plan | Where-Object { -not $done.ContainsKey((Key $_)) })
if ($pendF.Count -eq 0) {
    Say 'CAMPANHA COMPLETA. Proximo passo: regenerar tabelas e figuras (me chame).'
} else {
    Say ("ainda pendentes {0} runs; relancar com -Mode Apply retoma do ponto." -f $pendF.Count)
}
