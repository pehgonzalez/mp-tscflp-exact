#Requires -Version 5.1
<#
.SYNOPSIS
  Extensao da rodada 3, cinco blocos, com tela de acompanhamento ao vivo.

.DESCRIPTION
  Substitui a invocacao manual de

    scripts\run_campanha_revisao.ps1 -Mode Apply -Blocks "cert-dump,mip-fresco,fatorial-100"

  por um script dedicado, e acrescenta dois blocos novos que nao existem no
  script original (heavytail e fase-instr). A semantica de plano, retomada,
  journal, epoca, STOP, env var do certificado e build e a MESMA do original:
  o corpo do plano (Build-Plan), a chave de "ja feito" (Load-Done/Key) e o
  preambulo de build foram copiados sem alteracao, e os blocos novos so
  acrescentam itens ao plano, sem mexer nos que ja existiam. Para os tres
  blocos herdados, rodar este script ou o original com
  -Blocks "cert-dump,mip-fresco,fatorial-100" produz o mesmo conjunto de runs
  e os mesmos registros.

  Ordem de execucao: cert-dump, heavytail, fase-instr, mip-fresco,
  fatorial-100.

    bloco cert-dump      R2. As 12 fechadas de novo, cada uma com a engine que
                         a fechou no lote canonico, com MPTSCFL_CERT_DUMP
                         apontando para logs\certdump. Rotulo exact-cert3.
                         12 runs, 8 a 1800 s e as 4 que so fecham no orcamento
                         estendido a 14400 s. 20.0 h. O env var e setado so
                         neste bloco e limpo depois de cada run.
    bloco heavytail      R11. As duas caudas pesadas (PSC2-C4-100-5 e
                         PSC3-C4-100-5) em m2 exact-noldb, isto e, com o
                         cutoff e a linha do bound retirados. 2 runs, 1.0 h.
    bloco fase-instr     Trajetoria instrumentada [lag] da fase Lagrangiana.
                         As 100 instancias, m2, TL 270 s, dois runs cada:
                         exact (laco completo) e exact-oneeval (uma unica
                         avaliacao), no mesmo orcamento. 200 runs, 15.0 h. O
                         dado sai no STDOUT de cada run, preservado em
                         logs\campanha_revisao\stdout\.
    bloco mip-fresco     R7. O MIP compacto (m0, exact) nas 100 instancias a
                         1800 s, no build fresco. 100 runs, 50.0 h. O conjunto
                         50-10 ja rodou pos-epoca com esta chave exata e e
                         pulado sozinho, entao na pratica sao 75 runs.
    bloco fatorial-100   R8. O fatorial de 5 bracos do bloco fatorial-50-5
                         repetido no conjunto 100-10. 125 runs, 62.5 h. Os
                         bracos m1 exact e m2 exact coincidem com o
                         100-sem-corepoint e sao pulados se aquele bloco ja
                         rodou, o que deixa 75 runs.

  439 runs no plano dos cinco blocos, 148.5 h. Contra um results.csv que ja
  tenha o lote canonico pos-epoca, 75 sao pulados sozinhos.

  STDOUT POR RUN. O stdout de todo run e preservado em
  logs\campanha_revisao\stdout\<inst>_m<method>_<mode>_tl<TL>.out. E de la que
  saem as linhas [lag] da instrumentacao, o dado do bloco fase-instr. O
  arquivo e reescrito quando o mesmo run e refeito.

  CLOCKS. Em -Mode Apply sobe um job de telemetria que grava, a cada 60 s,
  uma linha em logs\campanha_revisao\clocks.csv com CurrentClockSpeed,
  MaxClockSpeed e a carga do _Total. E a documentacao de clock sustentado que
  o parecer pede (R15). O job morre junto com a campanha, inclusive por
  Ctrl-C. Onde o CIM nao responder, a linha sai com NA no campo e nada quebra.

  TELA. Em console interativo, o script reserva as primeiras linhas da janela
  e redesenha a cada 5 s um painel com uma linha por bloco, uma barra de
  progresso global e a ultima linha do .gurobi.log do run corrente. Fora de
  console interativo (saida redirecionada, RawUI ausente) cai sozinho para o
  modo texto simples, uma linha por evento, igual ao script original. Para
  forcar o modo texto, defina a variavel de ambiente MPTSCFL_NO_DASHBOARD.

  PARADA. Ctrl-C interrompe na hora (o run em andamento morre junto e a linha
  dele nao e gravada). O jeito limpo e criar o arquivo
  logs\campanha_revisao\STOP: o run corrente termina, e a campanha para antes
  do proximo. Relancar retoma do ponto exato nos dois casos.

.PARAMETER Mode
  Report  imprime o plano e o que ja esta feito, nao executa nada
  Apply   executa
.PARAMETER Threads
  Threads do Gurobi, repassado ao binario. Padrao 16, como no original.
.PARAMETER GurobiHome
  So e preciso quando a variavel de ambiente GUROBI_HOME nao existe.
  Exemplo: -GurobiHome C:\gurobi1302\win64
.EXAMPLE
  powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\run_extensao_rodada3.ps1 -Mode Report
  powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\run_extensao_rodada3.ps1 -Mode Apply
#>
[CmdletBinding()]
param(
    [ValidateSet('Report','Apply')]
    [string]$Mode = 'Report',
    [int]$Threads = 16,
    # So e preciso quando a variavel de ambiente GUROBI_HOME nao existe.
    # Exemplo: -GurobiHome C:\gurobi1302\win64
    [string]$GurobiHome = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version 2.0

# Os glifos do painel sao UTF-8; sem isto o console do Windows imprime lixo.
# (O proprio arquivo esta salvo em UTF-8 COM BOM: o PowerShell 5.1 le UTF-8
# sem BOM como ANSI e quebraria os mesmos glifos no codigo fonte.)
try { [Console]::OutputEncoding = [System.Text.Encoding]::UTF8 } catch { }

# A -Seed global do script original. Nenhum dos cinco blocos desta extensao usa
# outra seed (o bloco sementes, o unico com seed != 0, nao esta aqui), entao o
# valor fica fixo em 0 e a chave de "ja feito" bate com a do original.
$Seed = 0

# O binario grava logs\results.csv relativo ao diretorio corrente, entao o
# script se muda para a raiz do repositorio, a pasta pai de scripts\.
if ($PSScriptRoot) { Set-Location (Split-Path -Parent $PSScriptRoot) }
$Root = (Get-Location).Path

$Exe      = Join-Path $Root 'build\Release\mptscfl.exe'
$TestExe  = Join-Path $Root 'build\Release\test_core.exe'
$DataDir  = 'data\instances'
$Results  = Join-Path $Root 'logs\results.csv'
$LogDir   = Join-Path $Root 'logs'
$StateDir = Join-Path $Root 'logs\campanha_revisao'
$EpochF   = Join-Path $StateDir 'epoch.txt'
$StopF    = Join-Path $StateDir 'STOP'
$Journal  = Join-Path $StateDir 'journal.csv'
# Capturas do console do run corrente, usadas so quando o painel esta ligado
# (com o painel na tela o filho nao pode escrever direto no console). O .out e
# movido para stdout\ assim que o run termina; o .err so sobrevive nas falhas.
$OutF      = Join-Path $StateDir 'run_atual.out'
$ErrF      = Join-Path $StateDir 'run_atual.err'
# Stdout preservado por run (as linhas [lag] da instrumentacao saem daqui).
$StdoutDir = Join-Path $StateDir 'stdout'
# Telemetria de clock/carga, uma linha por minuto enquanto a campanha roda.
$ClocksF   = Join-Path $StateDir 'clocks.csv'

# Os cinco blocos desta extensao, na ordem em que Build-Plan os produz, que e a
# ordem de execucao.
$Blocos = @('cert-dump','heavytail','fase-instr','mip-fresco','fatorial-100')

function Say([string]$m) { Write-Host ("[{0}] {1}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $m) }

# ---------------------------------------------------------------- o plano ----
# Tudo desta secao e copia literal de scripts\run_campanha_revisao.ps1. O plano
# completo e construido igual e so depois filtrado pelos tres blocos, para nao
# haver chance de a ordem ou os campos divergirem do original.

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
    # ------------------------------------------ blocos novos desta extensao ----
    # Os dois blocos abaixo NAO existem no run_campanha_revisao.ps1. Entram aqui,
    # entre cert-dump e mip-fresco, porque a ordem de execucao pedida e
    # cert-dump, heavytail, fase-instr, mip-fresco, fatorial-100. Como so
    # acrescentam itens, os blocos herdados continuam com o mesmo conteudo e a
    # mesma ordem relativa do original.
    # heavytail (R11): as duas caudas pesadas com o cutoff e a linha do bound
    # retirados, que e exatamente o que o rotulo exact-noldb faz.
    & $add 'heavytail' 'PSC2-C4-100-5' 2 'exact-noldb' 1800 0
    & $add 'heavytail' 'PSC3-C4-100-5' 2 'exact-noldb' 1800 0
    # fase-instr: trajetoria [lag] da fase Lagrangiana nas 100 instancias, dois
    # runs por instancia no MESMO orcamento de 270 s, o laco completo (exact) e
    # a avaliacao unica (exact-oneeval), nesta ordem. O TL de 270 s e o que
    # separa a chave destes runs da dos blocos de 1800 s.
    foreach ($sfx in $Sets) {
        foreach ($i in (All-Instances $sfx)) {
            & $add 'fase-instr' $i 2 'exact'         270 0
            & $add 'fase-instr' $i 2 'exact-oneeval' 270 0
        }
    }
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

# ------------------------------------------------------------ o painel ----
# Estado unico do painel. Todas as chaves nascem inicializadas por causa do
# Set-StrictMode.

$script:D = @{
    On       = $false   # painel ligado (host interativo e -Mode Apply)
    W        = 0        # largura conhecida da janela
    H        = 0        # altura conhecida da janela
    Rows     = 11       # linhas reservadas no topo (2 separadores + 5 blocos +
                        # barra + cauda do log + mensagem + 1 em branco)
    T0       = (Get-Date)
    Total    = 0        # runs do plano (os tres blocos)
    Done     = 0        # runs efetivos concluidos (inclui os ja no results.csv)
    Ratio    = 1.0
    RemTL    = 0.0      # soma dos TL pendentes a partir do run corrente
    BTotal   = @{}      # bloco -> runs no plano
    BDone    = @{}      # bloco -> runs concluidos
    BRemTL   = @{}      # bloco -> soma dos TL pendentes do bloco
    BT0      = @{}      # bloco -> inicio nesta sessao
    BDur     = @{}      # bloco -> duracao medida nesta sessao (segundos)
    CurBloco = ''
    CurInst  = ''
    CurStart = $null
    CurLog   = ''
    Msg      = ''
    MsgColor = 'DarkGray'
}

function Get-ConSize {
    # Largura e altura uteis da janela, ou $null se o host nao souber dizer.
    # [Console] e a fonte primaria (e o que SetCursorPosition enxerga); em hosts
    # que reportam 0 (terminal sem tamanho negociado) cai para o RawUI.
    $w = 0; $h = 0
    try { $w = [int][Console]::WindowWidth; $h = [int][Console]::WindowHeight } catch { $w = 0; $h = 0 }
    if ($w -le 0 -or $h -le 0) {
        try {
            $ws = $Host.UI.RawUI.WindowSize
            if ($null -ne $ws) { $w = [int]$ws.Width; $h = [int]$ws.Height }
        } catch { }
    }
    if ($w -le 0 -or $h -le 0) { return $null }
    return @{ W = $w; H = $h }
}

function Test-Interactive {
    # Detectado UMA vez, no inicio. Qualquer duvida cai para o modo texto.
    try {
        if ($env:MPTSCFL_NO_DASHBOARD) { return $false }
        if ([Console]::IsOutputRedirected) { return $false }
        if ([Console]::IsErrorRedirected) { return $false }
        if ($null -eq $Host -or $null -eq $Host.UI -or $null -eq $Host.UI.RawUI) { return $false }
        $ws = $Host.UI.RawUI.WindowSize
        if ($null -eq $ws) { return $false }
        $sz = Get-ConSize
        if ($null -eq $sz) { return $false }
        if ($sz.W -lt 60 -or $sz.H -lt 16) { return $false }
        $null = [Console]::CursorTop
        return $true
    } catch { return $false }
}

function Fmt-HMS([double]$s) {
    if ([double]::IsNaN($s) -or [double]::IsInfinity($s) -or $s -lt 0) { $s = 0 }
    if ($s -gt 3599999) { $s = 3599999 }
    $t = [TimeSpan]::FromSeconds([Math]::Round($s))
    return ('{0}:{1:00}:{2:00}' -f [int][Math]::Floor($t.TotalHours), $t.Minutes, $t.Seconds)
}

function Fmt-HHMMSS([double]$s) {
    if ([double]::IsNaN($s) -or [double]::IsInfinity($s) -or $s -lt 0) { $s = 0 }
    if ($s -gt 3599999) { $s = 3599999 }
    $t = [TimeSpan]::FromSeconds([Math]::Round($s))
    return ('{0:00}:{1:00}:{2:00}' -f [int][Math]::Floor($t.TotalHours), $t.Minutes, $t.Seconds)
}

function Find-RunLog([string]$inst, $md, [string]$mode, $sd, [datetime]$since) {
    # logs\<inst>_m<method>_<mode>_s<seed>_<timestamp>.gurobi.log, o nome que o
    # main.cpp monta. Pega o mais recente que nao seja anterior ao run.
    try {
        $pat = '{0}_m{1}_{2}_s{3}_*.gurobi.log' -f $inst, $md, $mode, $sd
        $f = @(Get-ChildItem -Path $LogDir -Filter $pat -File -ErrorAction SilentlyContinue |
               Where-Object { $_.LastWriteTime -ge $since.AddSeconds(-10) } |
               Sort-Object LastWriteTime -Descending)
        if ($f.Count -gt 0) { return $f[0].FullName }
    } catch { }
    return ''
}

function Get-LogTail([string]$path, [int]$max) {
    # Ultima linha nao vazia. Abre compartilhado: o Gurobi mantem o arquivo
    # aberto para escrita enquanto resolve.
    if ([string]::IsNullOrEmpty($path)) { return '' }
    if ($max -lt 1) { return '' }
    try {
        if (-not (Test-Path -LiteralPath $path)) { return '' }
        $share = [System.IO.FileShare]::ReadWrite -bor [System.IO.FileShare]::Delete
        $fs = New-Object System.IO.FileStream($path, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, $share)
        try {
            $len  = $fs.Length
            if ($len -le 0) { return '' }
            $take = [int][Math]::Min([long]8192, $len)
            [void]$fs.Seek($len - $take, [System.IO.SeekOrigin]::Begin)
            $buf = New-Object byte[] $take
            $n   = $fs.Read($buf, 0, $take)
            $txt = [System.Text.Encoding]::UTF8.GetString($buf, 0, $n)
        } finally { $fs.Dispose() }
        $lines = $txt -split "`r?`n"
        for ($i = $lines.Count - 1; $i -ge 0; $i--) {
            $l = $lines[$i].Trim()
            if ($l -ne '') {
                if ($l.Length -gt $max) { return $l.Substring(0, $max) }
                return $l
            }
        }
    } catch { return '' }
    return ''
}

function Get-StdoutPath($r) {
    # logs\campanha_revisao\stdout\<inst>_m<method>_<mode>_tl<TL>.out, com os
    # caracteres que o Windows nao aceita em nome de arquivo trocados por '_'.
    # (Os rotulos usados aqui ja sao seguros; a limpeza e cinto de seguranca
    # para rotulos futuros.)
    $n = '{0}_m{1}_{2}_tl{3}.out' -f $r.Inst, $r.Method, $r.Md, $r.TL
    $n = $n -replace '[\\/:\*\?"<>\|]', '_'
    $n = $n -replace '[^A-Za-z0-9._\-]', '_'
    if ($n.Length -gt 120) { $n = $n.Substring(0, 120) }
    return (Join-Path $StdoutDir $n)
}

function New-Seg([string]$t, [string]$c) { return [pscustomobject]@{ T = $t; C = $c } }

function Write-Row([int]$row, $segs) {
    # Pinta uma linha inteira na posicao dada e completa com espacos ate a
    # penultima coluna (escrever na ultima coluna rola a tela).
    $w = $script:D.W
    if ($w -lt 20) { return }
    try { [Console]::SetCursorPosition(0, $row) } catch { return }
    $used = 0
    foreach ($s in $segs) {
        if ($used -ge ($w - 1)) { break }
        $t = [string]$s.T
        if ($t.Length -gt ($w - 1 - $used)) { $t = $t.Substring(0, $w - 1 - $used) }
        if ($t.Length -gt 0) {
            try { Write-Host -NoNewline $t -ForegroundColor $s.C } catch { Write-Host -NoNewline $t }
            $used += $t.Length
        }
    }
    if ($used -lt ($w - 1)) { Write-Host -NoNewline (' ' * (($w - 1) - $used)) }
}

function Pad-LR([string]$l, [string]$r, [int]$w) {
    # Junta esquerda e direita numa linha de largura $w, truncando a esquerda.
    if ($w -lt 4) { return '' }
    if ([string]::IsNullOrEmpty($r)) {
        if ($l.Length -gt $w) { return $l.Substring(0, $w) }
        return $l
    }
    $gap = $w - $l.Length - $r.Length
    if ($gap -ge 2) { return ($l + (' ' * $gap) + $r) }
    # Janela estreita: a coluna da direita e sacrificada, o lado esquerdo (o
    # progresso do bloco) e o que interessa.
    if ($l.Length -gt $w) { return $l.Substring(0, $w) }
    return $l
}

function Start-Dash {
    if (-not $script:D.On) { return }
    $sz = Get-ConSize
    if ($null -eq $sz) { $script:D.On = $false; return }
    $script:D.W = $sz.W
    $script:D.H = $sz.H
    try { [Console]::CursorVisible = $false } catch { }
    try { [Console]::Clear() } catch { $script:D.On = $false; return }
    Show-Dash
}

function Stop-Dash {
    if (-not $script:D.On) { return }
    Show-Dash
    try {
        [Console]::SetCursorPosition(0, [Math]::Min($script:D.Rows, $script:D.H - 1))
        [Console]::CursorVisible = $true
    } catch { }
    Write-Host ''
}

function Show-Dash {
    if (-not $script:D.On) { return }
    $sz = Get-ConSize
    if ($null -eq $sz) { $script:D.On = $false; return }
    $w = $sz.W; $h = $sz.H
    # Janela redimensionada: limpa e repinta tudo com as larguras novas.
    if ($w -ne $script:D.W -or $h -ne $script:D.H) {
        $script:D.W = $w; $script:D.H = $h
        try { [Console]::Clear() } catch { }
    }
    # Estreita ou baixa demais para o painel: nao pinta nada nesta rodada.
    if ($w -lt 50 -or $h -lt ($script:D.Rows + 1)) { return }

    $now    = Get-Date
    $sessEl = ($now - $script:D.T0).TotalSeconds
    $curEl  = 0.0
    if ($null -ne $script:D.CurStart) { $curEl = ($now - $script:D.CurStart).TotalSeconds }

    $row = 0
    Write-Row $row @((New-Seg ('-' * ($w - 1)) 'DarkGray')); $row++

    foreach ($b in $Blocos) {
        $tot = [int]$script:D.BTotal[$b]
        $dn  = [int]$script:D.BDone[$b]
        if ($b -eq $script:D.CurBloco) {
            $rem = ([double]$script:D.BRemTL[$b] * $script:D.Ratio) - $curEl
            if ($rem -lt 0) { $rem = 0 }
            $left  = ('  > ' + $b).PadRight(24) + ('{0} de {1} runs, prevista em ~{2}' -f $dn, $tot, (Fmt-HMS $rem))
            $ini   = $script:D.BT0[$b]
            $eb    = 0.0
            if ($null -ne $ini) { $eb = ($now - $ini).TotalSeconds }
            $right = 'inicio {0:HH:mm:ss}, decorrido {1}' -f $ini, (Fmt-HMS $eb)
            $cl = 'Yellow'; $cr = 'Cyan'
        } elseif ($dn -ge $tot -and $tot -gt 0) {
            $left = ('  ' + [string][char]0x2713 + ' ' + $b).PadRight(24) + ('{0} runs' -f $tot)
            if ($script:D.BDur.ContainsKey($b) -and $null -ne $script:D.BDur[$b]) {
                $right = 'concluida em {0}' -f (Fmt-HMS ([double]$script:D.BDur[$b]))
            } else {
                $right = 'concluida antes desta sessao'
            }
            $cl = 'Green'; $cr = 'DarkGray'
        } else {
            $left  = ('  ' + [string][char]0x00B7 + ' ' + $b).PadRight(24) + ('{0} runs' -f $tot)
            $right = ''
            if ($dn -gt 0) { $right = ('{0} de {1} ja no results.csv' -f $dn, $tot) }
            $cl = 'DarkGray'; $cr = 'DarkGray'
        }
        $line = Pad-LR $left $right ($w - 1)
        # Duas cores: o corpo do bloco e a coluna da direita.
        if ($right -ne '' -and $line.EndsWith($right)) {
            $head = $line.Substring(0, $line.Length - $right.Length)
            Write-Row $row @((New-Seg $head $cl), (New-Seg $right $cr))
        } else {
            Write-Row $row @((New-Seg $line $cl))
        }
        $row++
    }

    Write-Row $row @((New-Seg ('-' * ($w - 1)) 'DarkGray')); $row++

    # Barra global. ~70 colunas, encolhida se a janela nao couber.
    $frac = 0.0
    if ($script:D.Total -gt 0) { $frac = [double]$script:D.Done / [double]$script:D.Total }
    if ($frac -lt 0) { $frac = 0.0 }
    if ($frac -gt 1) { $frac = 1.0 }
    $tail = '  {0,5:N1}%  decorrido {1}  falta ~{2}'
    $rem  = ([double]$script:D.RemTL * $script:D.Ratio) - $curEl
    if ($rem -lt 0) { $rem = 0 }
    $tailS = $tail -f ($frac * 100), (Fmt-HHMMSS $sessEl), (Fmt-HHMMSS $rem)
    $barW = $w - 2 - $tailS.Length
    if ($barW -gt 70) { $barW = 70 }
    if ($barW -lt 10) { $barW = 10 }
    $fill = [int][Math]::Floor($frac * $barW)
    if ($fill -gt $barW) { $fill = $barW }
    $segs = @()
    $segs += (New-Seg ' ' 'Gray')
    if ($fill -gt 0)          { $segs += (New-Seg ([string][char]0x2588 * $fill) 'Green') }
    if (($barW - $fill) -gt 0){ $segs += (New-Seg ([string][char]0x2591 * ($barW - $fill)) 'DarkGray') }
    $segs += (New-Seg $tailS 'White')
    Write-Row $row $segs; $row++

    # Cauda do log do run corrente.
    $t = ''
    if ($script:D.CurLog -ne '') { $t = Get-LogTail $script:D.CurLog ($w - 5) }
    if ($t -eq '' -and $script:D.CurInst -ne '') { $t = ('{0} ... aguardando a primeira linha do log' -f $script:D.CurInst) }
    Write-Row $row @((New-Seg '  ' 'DarkGray'), (New-Seg $t 'DarkGray')); $row++

    # Ultimo aviso, ou a linha de ajuda.
    $m = $script:D.Msg
    if ($m -eq '') { $m = ('STOP limpo: crie ' + $StopF) }
    Write-Row $row @((New-Seg '  ' 'DarkGray'), (New-Seg $m $script:D.MsgColor)); $row++

    Write-Row $row @((New-Seg '' 'Gray'))
    try { [Console]::SetCursorPosition(0, [Math]::Min($script:D.Rows, $h - 1)) } catch { }
}

function Set-DashMsg([string]$m, [string]$color) {
    $script:D.Msg = $m
    $script:D.MsgColor = $color
}

# ----------------------------------------------------------------- inicio ----

Write-Host ('=' * 78)
Write-Host ("Extensao da rodada 3 ({0}), modo {1}" -f ($Blocos -join ', '), $Mode)
Write-Host ('=' * 78)

if (-not (Test-Path -LiteralPath (Join-Path $Root $DataDir))) {
    throw "data\instances nao encontrada em $Root. Rode da raiz do repositorio."
}
$nInst = @(Get-ChildItem -Path (Join-Path $Root $DataDir) -Filter 'PSC*.txt' -File).Count
if ($nInst -ne 100) { throw "data\instances tem $nInst instancias, esperava 100." }

$plan = Build-Plan
$plan = @($plan | Where-Object { $Blocos -contains $_.Bloco })

if (-not (Test-Path -LiteralPath $StateDir)) { New-Item -ItemType Directory -Path $StateDir -Force | Out-Null }

# Epoca. Definida uma unica vez, DEPOIS do primeiro build fresco. Nunca criada
# aqui se ja existe: e a mesma epoca da campanha inteira.
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
    if (-not (Test-Path -LiteralPath $StdoutDir)) { New-Item -ItemType Directory -Path $StdoutDir -Force | Out-Null }
}
if ($epoch -eq '') { $epoch = '99999999-999999' } # Report antes do 1o Apply: nada feito

$done = Load-Done $epoch
$pend = @($plan | Where-Object { -not $done.ContainsKey((Key $_)) })
$feitos = $plan.Count - $pend.Count

$somaTL = 0.0
foreach ($r in $pend) { $somaTL += [double]$r.TL }

Write-Host ''
Say ("plano total {0} runs | concluidos {1} | pendentes {2}" -f $plan.Count, $feitos, $pend.Count)
foreach ($b in $Blocos) {
    $bg = @($plan | Where-Object { $_.Bloco -eq $b })
    $pg = @($pend | Where-Object { $_.Bloco -eq $b })
    $h = 0.0; foreach ($r in $pg) { $h += [double]$r.TL }
    Write-Host ("  {0,-20} {1,3} runs, pendentes {2,3}, ~{3,6:N1} h" -f $b, $bg.Count, $pg.Count, ($h/3600))
}
Say ("orcamento pendente somado: {0:N1} h (~{1:N1} dias)" -f ($somaTL/3600), ($somaTL/86400))

if ($Mode -eq 'Report') {
    Write-Host ''
    Write-Host 'Nada foi executado. Para rodar:' -ForegroundColor Green
    Write-Host '  powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\run_extensao_rodada3.ps1 -Mode Apply' -ForegroundColor Green
    Write-Host ''
    Write-Host 'Durante o Apply:' -ForegroundColor Green
    Write-Host ('  parada limpa, crie o arquivo ' + $StopF) -ForegroundColor Green
    Write-Host '  (o run em andamento termina e a campanha para antes do proximo);' -ForegroundColor Green
    Write-Host '  Ctrl-C tambem para, mas mata o run em andamento sem gravar a linha dele.' -ForegroundColor Green
    Write-Host '  Relancar com -Mode Apply retoma do ponto exato nos dois casos.' -ForegroundColor Green
    Write-Host ('  stdout de cada run preservado em ' + $StdoutDir) -ForegroundColor Green
    Write-Host ('  clock/carga amostrados a cada 60 s em ' + $ClocksF) -ForegroundColor Green
    return
}

# ------------------------------------------------------------ o laco real ----

$ratio = 1.0      # razao media tempo real / TL, refinada a cada run
$nMed  = 0
$seq   = 0
$falhas = 0
$campT0 = Get-Date

# Telemetria de clock (R15). Um job em segundo plano grava uma linha por minuto
# enquanto a campanha roda; e o registro de que o clock sustentado nao caiu no
# meio do lote. O job e um processo PowerShell separado, entao recebe o caminho
# ABSOLUTO por argumento (jobs nao herdam o diretorio corrente) e nao enxerga
# nenhuma variavel do script. Tudo la dentro e try/catch: numa maquina em que o
# CIM nao responda a linha sai com NA e a campanha segue.
$ClocksJob = $null
try {
    $ClocksJob = Start-Job -Name 'mptscfl-clocks' -ArgumentList $ClocksF -ScriptBlock {
        param($path)
        if (-not (Test-Path -LiteralPath $path)) {
            try {
                Set-Content -LiteralPath $path -Encoding UTF8 `
                    -Value 'datetime,current_mhz,max_mhz,cpu_pct'
            } catch { }
        }
        while ($true) {
            $cur = 'NA'; $max = 'NA'; $pct = 'NA'
            try {
                $cpu = @(Get-CimInstance -ClassName Win32_Processor -ErrorAction Stop)
                if ($cpu.Count -gt 0) {
                    if ($null -ne $cpu[0].CurrentClockSpeed) { $cur = [string]$cpu[0].CurrentClockSpeed }
                    if ($null -ne $cpu[0].MaxClockSpeed)     { $max = [string]$cpu[0].MaxClockSpeed }
                    if ($null -ne $cpu[0].LoadPercentage)    { $pct = [string]$cpu[0].LoadPercentage }
                }
            } catch { }
            # O contador formatado do _Total e mais fiel que LoadPercentage e,
            # ao contrario de Get-Counter, nao depende do nome localizado do
            # contador ("\Processador(_Total)\..." em Windows pt-BR).
            try {
                $tot = @(Get-CimInstance -ClassName Win32_PerfFormattedData_PerfOS_Processor `
                            -Filter "Name='_Total'" -ErrorAction Stop)
                if ($tot.Count -gt 0 -and $null -ne $tot[0].PercentProcessorTime) {
                    $pct = [string]$tot[0].PercentProcessorTime
                }
            } catch { }
            try {
                Add-Content -LiteralPath $path -Encoding UTF8 -Value `
                    ('{0:yyyy-MM-ddTHH:mm:ss},{1},{2},{3}' -f (Get-Date), $cur, $max, $pct)
            } catch { }
            Start-Sleep -Seconds 60
        }
    }
    Say ("telemetria de clock ligada (job {0}), gravando em {1}" -f $ClocksJob.Id, $ClocksF)
} catch {
    $ClocksJob = $null
    Say ("AVISO: nao consegui subir o job de telemetria de clock ({0}); seguindo sem clocks.csv" -f $_.Exception.Message)
}

# Painel: detectado UMA vez, aqui.
$script:D.On    = (Test-Interactive)
$script:D.T0    = $campT0
$script:D.Total = $plan.Count
$script:D.Done  = $feitos
$script:D.Ratio = $ratio
$script:D.RemTL = $somaTL
foreach ($b in $Blocos) {
    $bg = @($plan | Where-Object { $_.Bloco -eq $b })
    $bp = @($pend | Where-Object { $_.Bloco -eq $b })
    $script:D.BTotal[$b] = $bg.Count
    $script:D.BDone[$b]  = $bg.Count - $bp.Count
    $s = 0.0; foreach ($r in $bp) { $s += [double]$r.TL }
    $script:D.BRemTL[$b] = $s
    $script:D.BT0[$b]    = $null
    $script:D.BDur[$b]   = $null
}
if (-not $script:D.On) {
    Say 'host nao interativo, seguindo em modo texto simples (uma linha por evento).'
}
Start-Dash

$avisos = New-Object System.Collections.ArrayList

# O finally e o que devolve o cursor e a tela ao normal quando o laco termina
# por Ctrl-C: sem ele o cursor ficaria invisivel no console do usuario.
try {
for ($idx = 0; $idx -lt $pend.Count; $idx++) {
    $r = $pend[$idx]

    if (Test-Path -LiteralPath $StopF) {
        if ($script:D.On) { Set-DashMsg 'STOP encontrado, parando com seguranca. Relancar retoma daqui.' 'Yellow'; Show-Dash }
        else { Say 'STOP encontrado, parando com seguranca. Relancar retoma daqui.' }
        break
    }
    # Alguem pode ter completado este run em outra janela; reconfere barato.
    $done2 = Load-Done $epoch
    if ($done2.ContainsKey((Key $r))) {
        if ($script:D.On) {
            $script:D.Done++
            $script:D.BDone[$r.Bloco] = [int]$script:D.BDone[$r.Bloco] + 1
            $script:D.BRemTL[$r.Bloco] = [double]$script:D.BRemTL[$r.Bloco] - [double]$r.TL
            $script:D.RemTL = [double]$script:D.RemTL - [double]$r.TL
            Set-DashMsg ("pulado (ja no results.csv): {0} {1} m{2}" -f $r.Inst, $r.Md, $r.Method) 'DarkGray'
            Show-Dash
        } else {
            Say ("pulado (ja no results.csv): {0} {1} m{2}" -f $r.Inst, $r.Md, $r.Method)
        }
        continue
    }

    $restante = 0.0
    for ($j = $idx; $j -lt $pend.Count; $j++) { $restante += [double]$pend[$j].TL }
    $etaS = $restante * $ratio
    $eta  = (Get-Date).AddSeconds($etaS)
    $entrada = Get-Date

    if ($script:D.On) {
        # Estado do painel para este run.
        $script:D.RemTL = $restante
        $rb = 0.0
        for ($j = $idx; $j -lt $pend.Count; $j++) {
            if ($pend[$j].Bloco -eq $r.Bloco) { $rb += [double]$pend[$j].TL }
        }
        $script:D.BRemTL[$r.Bloco] = $rb
        if ($null -eq $script:D.BT0[$r.Bloco]) { $script:D.BT0[$r.Bloco] = $entrada }
        $script:D.CurBloco = $r.Bloco
        $script:D.CurInst  = $r.Inst
        $script:D.CurStart = $entrada
        $script:D.CurLog   = ''
        $script:D.Ratio    = $ratio
        # Mensagem do run anterior nao vale mais para este.
        Set-DashMsg '' 'DarkGray'
        Show-Dash
    } else {
        Write-Host ''
        Write-Host ('-' * 78)
        Say ("[{0}] run {1}/{2} da fila | {3} | m{4} {5} | TL {6} s | seed {7}" -f `
            $r.Bloco, ($idx+1), $pend.Count, $r.Inst, $r.Method, $r.Md, $r.TL, $r.Seed)
        Say ("entrada {0:HH:mm:ss} | restante ~{1:N1} h | termino estimado da campanha {2:ddd dd/MM HH:mm}" -f `
            $entrada, ($etaS/3600), $eta)
    }

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
        if ($script:D.On) { Set-DashMsg ("MPTSCFL_CERT_DUMP={0} (so neste bloco)" -f $certDir) 'DarkGray' }
        else { Say ("MPTSCFL_CERT_DUMP={0} (so neste bloco)" -f $certDir) }
    }

    $instFile = Join-Path $DataDir ($r.Inst + '.txt')
    $rc = 0
    # Stdout preservado por run. Com o painel ligado o filho ja escreve num
    # arquivo, e so mover no fim; no modo texto a saida continua aparecendo no
    # console linha a linha (Write-Host no meio do pipe) e uma copia vai para o
    # mesmo destino. E daqui que sai a trajetoria [lag] do bloco fase-instr.
    $runOut = Get-StdoutPath $r
    if (-not (Test-Path -LiteralPath $StdoutDir)) {
        New-Item -ItemType Directory -Path $StdoutDir -Force | Out-Null
    }
    try {
        if ($script:D.On) {
            # Com o painel na tela o filho nao pode escrever no console, entao a
            # invocacao vira assincrona com redirecionamento: o laco de espera
            # (Wait de 5 s) e o que da o redesenho durante o run. O env var, os
            # argumentos e o exit code sao os mesmos do modo texto.
            $argl = @()
            foreach ($a in @($instFile, $r.Method, $r.TL, $r.Md, $r.Seed, $Threads, $r.Pap)) {
                $q = [string]$a
                if ($q -match '[\s"]') { $q = '"' + ($q -replace '"', '\"') + '"' }
                $argl += $q
            }
            $p = Start-Process -FilePath $Exe -ArgumentList ($argl -join ' ') -NoNewWindow -PassThru `
                    -RedirectStandardOutput $OutF -RedirectStandardError $ErrF
            while (-not $p.WaitForExit(5000)) {
                if ($script:D.CurLog -eq '') {
                    $script:D.CurLog = Find-RunLog $r.Inst $r.Method $r.Md $r.Seed $entrada
                    if ($script:D.CurLog -eq '' -and (Test-Path -LiteralPath $OutF)) { $script:D.CurLog = $OutF }
                }
                Show-Dash
            }
            $p.WaitForExit()
            $rc = $p.ExitCode
            if ($null -eq $rc) { $rc = -1 }
        } else {
            # Igual ao original, com o stdout passando por um pipe que reimprime
            # cada linha no console e grava a copia. O stderr continua indo
            # direto para o console, e $LASTEXITCODE continua sendo o do
            # binario mesmo com o pipeline no meio.
            & $Exe $instFile $r.Method $r.TL $r.Md $r.Seed $Threads $r.Pap |
                ForEach-Object { Write-Host $_; $_ } |
                Set-Content -LiteralPath $runOut -Encoding UTF8
            $rc = $LASTEXITCODE
        }
    } finally {
        if ($certOn -and (Test-Path Env:MPTSCFL_CERT_DUMP)) { Remove-Item Env:MPTSCFL_CERT_DUMP }
    }
    $dur = ((Get-Date) - $entrada).TotalSeconds

    # Com o painel, o stdout do run acabou de ser fechado; vira arquivo do run.
    if ($script:D.On) {
        try {
            if (Test-Path -LiteralPath $OutF) { Move-Item -LiteralPath $OutF -Destination $runOut -Force }
        } catch { }
    }

    $seq++
    Add-Content -LiteralPath $Journal -Encoding UTF8 -Value `
        ('{0},{1:yyyy-MM-ddTHH:mm:ss},{2},{3},{4},{5},{6},{7:N0},{8},{9}' -f `
        $seq, $entrada, $r.Bloco, $r.Inst, $r.Method, $r.Md, $r.TL, $dur, $rc,
        ('executado seed=' + $r.Seed))
    # A seed vai dentro da coluna acao de proposito: acrescentar uma coluna
    # nova quebraria os journal.csv ja gravados com o cabecalho antigo.

    if ($rc -ne 0) {
        $falhas++
        $av = ("AVISO: exit {0} em {1} m{2} {3}; seguindo (o run pode ser refeito relancando)" -f $rc, $r.Inst, $r.Method, $r.Md)
        [void]$avisos.Add($av)
        if ($script:D.On) {
            # Sem console para o filho, guarda a captura da falha para depois.
            # O stdout ja virou $runOut logo acima; aqui se junta o stderr.
            try {
                $fal = Join-Path $StateDir ('falha_{0}_m{1}_{2}_{3:yyyyMMdd-HHmmss}.log' -f $r.Inst, $r.Method, $r.Md, $entrada)
                $txt = ''
                if (Test-Path -LiteralPath $runOut) { $txt += (Get-Content -LiteralPath $runOut -Raw) }
                if (Test-Path -LiteralPath $ErrF)   { $txt += "`n--- stderr ---`n" + (Get-Content -LiteralPath $ErrF -Raw) }
                Set-Content -LiteralPath $fal -Value $txt -Encoding UTF8
            } catch { }
            Set-DashMsg $av 'Red'
        } else {
            Say $av
        }
    }
    # Refina a razao real/TL com media movel (runs fechados cedo puxam para baixo).
    $obs = $dur / [double]$r.TL
    if ($obs -gt 0.01) {
        $nMed++
        $ratio = (($ratio * ($nMed - 1)) + [Math]::Min($obs, 1.2)) / $nMed
        if ($ratio -lt 0.25) { $ratio = 0.25 }
    }

    if ($script:D.On) {
        $script:D.Ratio = $ratio
        $script:D.Done++
        $script:D.BDone[$r.Bloco] = [int]$script:D.BDone[$r.Bloco] + 1
        $script:D.BRemTL[$r.Bloco] = [double]$script:D.BRemTL[$r.Bloco] - [double]$r.TL
        $script:D.RemTL = [double]$script:D.RemTL - [double]$r.TL
        if ([int]$script:D.BDone[$r.Bloco] -ge [int]$script:D.BTotal[$r.Bloco]) {
            $ini = $script:D.BT0[$r.Bloco]
            if ($null -ne $ini) { $script:D.BDur[$r.Bloco] = ((Get-Date) - $ini).TotalSeconds }
            $script:D.CurBloco = ''
        }
        $script:D.CurInst  = ''
        $script:D.CurStart = $null
        $script:D.CurLog   = ''
        Show-Dash
    }
}
} finally {
    if ($script:D.On) {
        $script:D.CurBloco = ''
        $script:D.CurInst  = ''
        $script:D.CurStart = $null
        Stop-Dash
    }
    # A telemetria morre junto com a campanha, inclusive por Ctrl-C.
    if ($null -ne $ClocksJob) {
        try { Stop-Job   -Job $ClocksJob -ErrorAction SilentlyContinue } catch { }
        try { Remove-Job -Job $ClocksJob -Force -ErrorAction SilentlyContinue } catch { }
        $ClocksJob = $null
    }
}

Write-Host ''
Write-Host ('=' * 78)
$dTot = ((Get-Date) - $campT0)
Say ("sessao encerrada | runs executados {0} | falhas {1} | duracao {2:N1} h" -f $seq, $falhas, $dTot.TotalHours)
# No modo texto os avisos ja sairam na hora, como no original; com o painel na
# tela eles ficaram so na linha de mensagem, entao vao todos aqui no fim.
if ($script:D.On) { foreach ($a in $avisos) { Say $a } }
$done = Load-Done $epoch
$pendF = @($plan | Where-Object { -not $done.ContainsKey((Key $_)) })
Say ("stdout dos runs em {0} | telemetria em {1}" -f $StdoutDir, $ClocksF)
if ($pendF.Count -eq 0) {
    Say ("EXTENSAO DA RODADA 3 COMPLETA ({0})." -f ($Blocos -join ', '))
} else {
    Say ("ainda pendentes {0} runs; relancar com -Mode Apply retoma do ponto." -f $pendF.Count)
}
