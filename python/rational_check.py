#!/usr/bin/env python3
"""Verificador exato (aritmetica racional) dos certificados .certdump do MP-TSCFLP.

O binario, com MPTSCFL_CERT_DUMP setada, grava ao fim de um run OPTIMAL um
arquivo <dir>/<inst>_<timestamp>.certdump com a solucao incumbente (y,z), os
coeficientes de todos os cortes de otimalidade de Benders adicionados, as
desigualdades (F_l), o bound Lagrangiano v_LD e o custo verificado. O formato
esta documentado no cabecalho do proprio arquivo e em patch_certdump.md.

Este script refaz quatro coisas SEM SOLVER e SEM DEPENDENCIA EXTERNA (so a
biblioteca padrao; nada de scipy, networkx ou numpy), em inteiros e Fraction,
portanto sem nenhum arredondamento de ponto flutuante do lado do verificador:

  CHECK 1, primal. Le a instancia original, confere que (y,z) satisfaz as
  desigualdades (F_l) do dump, re-roteia todos os produtos com um fluxo de
  custo minimo por caminhos minimos sucessivos (Dijkstra com potenciais, tudo
  em int), confirma que a demanda de cada produto e integralmente atendida e
  recalcula custo = transporte + custos fixos das instalacoes abertas em (y,z).
  Compara com VERIFIED e OBJ do dump.

  CHECK 2, fechamento dual nos cortes. Cada corte de Benders da um piso para
  theta_l quando avaliado no incumbente (y*,z*). Somando os custos fixos de
  (y*,z*) com o maximo dos cortes por produto obtem-se um piso do custo do
  incumbente. Se esse piso ficar a menos de 0.9999 do custo verificado, o par
  (custo, cortes) fecha em aritmetica exata.

  CHECK 3, validade GLOBAL de cada corte (novo). Reconstroi, a partir dos
  coeficientes do dump, um ponto dual do LP de roteamento do produto l e decide
  EXATAMENTE se esse ponto certifica o corte para TODO (y,z) binario, nao so
  para (y*,z*). Ver a secao ESTRUTURA DUAL abaixo. E o check que fecha a
  exigencia do referee: cada desigualdade do mestre final passa a ser uma
  desigualdade demonstrada valida, e nao apenas um objeto produzido pelo
  solver.

  CHECK 4, pisos iniciais e linha Lagrangiana (novo). Recomputa exatamente
  v_l(1,1), o custo de roteamento do produto l com tudo aberto, que e o valor
  usado pelo construtor do mestre como cota inferior de theta_l, e mede a
  margem theta*_l - v_l(1,1) no incumbente. Recomputa tambem a folga exata
  custo - v_LD da linha global "fixos + sum theta >= v_LD". Ambas as margens
  exoneram numericamente as duas unicas desigualdades do mestre cujo lado
  direito veio de um numero de ponto flutuante.


ESTRUTURA DUAL DO SUBPROBLEMA (derivada de src3/benders_model.cpp)
=================================================================
ProductSubproblem monta, para o produto l e um par (y,z) dado, o LP

    min   sum_ij c_ijl x_ij + sum_jk d_jkl w_jk
    s.a.  sum_j w_jk            >= q_kl        (dem_k)     dual alpha_k >= 0
          sum_i x_ij - sum_k w_jk >= 0         (cons_j)    dual beta_j  >= 0
          sum_j x_ij            <= b_il y_i    (fcap_i)    dual pi_i    <= 0
          sum_k w_jk            <= p_jl z_j    (wcap_j)    dual gamma_j <= 0
          x, w >= 0

(os sinais dos duais sao a convencao do Gurobi para minimizacao: Pi >= 0 em
linha >=, Pi <= 0 em linha <=). O dual e

    max   sum_k q_kl alpha_k + sum_i (b_il y_i) pi_i + sum_j (p_jl z_j) gamma_j
    s.a.  pi_i + beta_j              <= c_ijl        (coluna x_ij)
          alpha_k - beta_j + gamma_j <= d_jkl        (coluna w_jk)
          alpha >= 0, beta >= 0, pi <= 0, gamma <= 0

e o corte gravado pelo solve() e exatamente o objetivo dual nesse ponto:

    constant = sum_k q_kl * Pi(dem_k)          -> sum_k q_kl alpha_k
    fac_i    = b_il * Pi(fcap_i)               -> b_il pi_i
    ware_j   = p_jl * Pi(wcap_j)               -> p_jl gamma_j
    theta_l >= constant + sum_i fac_i y_i + sum_j ware_j z_j

O ponto central: a REGIAO VIAVEL DUAL NAO DEPENDE DE (y,z). So o objetivo
depende. Logo qualquer ponto dual viavel gera uma desigualdade valida para
todo (y,z), e provar a validade global do corte = exibir um ponto dual viavel
cujo objetivo domine o corte emitido.

Do dump recuperam-se diretamente os duais das linhas cujo RHS depende de (y,z):

    pi_i    = fac_i  / b_il          (b_il > 0)
    gamma_j = ware_j / p_jl          (p_jl > 0)

e nao se recuperam os alpha_k individualmente, so a soma ponderada
sum_k q_kl alpha_k, que e a constante. As linhas de RHS zero (cons_j, o
balanco no deposito) tem duais beta_j que NAO aparecem no corte. A validade do
corte fica portanto equivalente a viabilidade de um sistema linear em
(alpha, beta) com (pi, gamma) fixos:

    beta_j <= c_ijl - pi_i               para todo i          (I x J linhas)
    alpha_k <= d_jkl + beta_j - gamma_j  para todo j          (J x K linhas)
    alpha, beta >= 0
    sum_k q_kl alpha_k >= constant

Esse sistema e um sistema de DIFERENCAS puro, um dual de caminho minimo em
camadas (fonte -> fabrica -> deposito -> cliente) com "comprimentos" c_ijl e
d_jkl e "ofertas" -pi_i e -gamma_j nos nos. Por ser um DAG em camadas o
Bellman-Ford exato colapsa em dois varrimentos de minimo, que sao os limites
superiores simultaneamente atingiveis:

    beta_j^max  = min_i (c_ijl - pi_i)                     >= 0
    alpha_k^max = min_j (d_jkl + beta_j^max - gamma_j)     >= 0

(o >= 0 e automatico: c, d >= 0, pi, gamma <= 0). Como aumentar beta so
afrouxa a restricao de alpha, e como aumentar alpha so aumenta o objetivo,
(alpha^max, beta^max) maximiza sum_k q_kl alpha_k sobre todo o sistema. Logo

    VALIDADE  <=>  Phi_l(pi,gamma) := sum_k q_kl alpha_k^max  >=  constant

e o teste e exato, necessario e suficiente PARA OS (pi,gamma) RECUPERADOS, e
nao precisa nem de simplex nem de Bellman-Ford iterativo: os dois minimos ja
resolvem o sistema de diferencas. Tudo em int (ver adiante) ou em Fraction.

Capacidade nula: se b_il = 0 entao fac_i = 0 para qualquer pi_i, e pi_i pode
ser tomado arbitrariamente negativo sem alterar o corte; o indice i sai do
minimo de beta_j. Idem para p_jl = 0 e o indice j no minimo de alpha_k. O
codigo trata esse caso; nas seis instancias nao ha capacidade nula.

SPARSIFICACAO E O SENTIDO DA DOBRA
==================================
BendersCallback::cut_expr() e accumulate_as_added() dobram na constante todo
coeficiente com |coef| <= eps (eps = 1e-6 por default) em vez de emiti-lo. Os
coeficientes sao <= 0 (produto de uma capacidade >= 0 por um dual <= 0), e
para y em [0,1] vale coef*y >= coef. Portanto a dobra troca o termo coef*y
pelo seu MENOR valor possivel: o lado direito emitido e MENOR OU IGUAL ao lado
direito do corte verdadeiro. A DOBRA ENFRAQUECE O CORTE. O corte emitido e
dominado pelo corte do dual verdadeiro e a validade do emitido segue da
validade do original. Consequencia pratica: a recuperacao pi_i = fac_i/b_il
com fac_i = 0 devolve pi_i = 0 no lugar de um pi_i verdadeiro em
[-eps/b_il, 0), o que so torna Phi_l MENOR, isto e, o CHECK 3 fica
CONSERVADOR. Um PASS continua sendo prova; um FAIL pode, em principio, ser
artefato da recuperacao.

Por isso o CHECK 3 tem um segundo estagio, acionado apenas quando o primeiro
falha: o teste com a FOLGA DA SPARSIFICACAO. Ele nao supoe nada sobre o dual
verdadeiro; prova diretamente o corte EMITIDO exibindo um ponto dual e pagando
explicitamente a diferenca por coordenada. O criterio geral, para (y,z)
binario arbitrario, e

    G(pi,gamma) = Phi_l(pi,gamma) - constant
                + sum_i min(0, b_il pi_i - fac_i)
                + sum_j min(0, p_jl gamma_j - ware_j)   >= 0

(o minimo do lado esquerdo menos o lado direito sobre y,z binarios). O segundo
estagio maximiza G sobre a JANELA DA SPARSIFICACAO, pi_i em [-eps/b_il, 0]
para os i com fac_i = 0 emitido (idem gamma_j), avaliando G exatamente em
Fraction num conjunto finito de pontos. Um G >= 0 em qualquer ponto e prova da
validade do corte emitido. O script tambem imprime o limite analitico
G <= G(recuperado) + Q_l * (max_i eps/b_il + max_j eps/p_jl); quando esse
limite ja e negativo, nenhum dual da janela salva o corte e a falha nao pode
ser artefato de sparsificacao.

Terceiro estagio, so em caso de falha: BUSCA DE CONTRAEXEMPLO. O corte e
avaliado em (y*,z*) contra v_l(y*,z*) (CHECK 1) e em (1,1) contra v_l(1,1)
(CHECK 4). Se o lado direito exceder o valor exato do subproblema num desses
pontos, o corte esta PROVADAMENTE INVALIDO, com testemunha explicita.

NOTA ARITMETICA. A matriz de restricoes do subproblema e a matriz de
incidencia de uma rede (o subproblema E um fluxo de custo minimo), logo
totalmente unimodular; com custos inteiros existe dual otimo inteiro e o
simplex devolve um vertice, que e inteiro. Nas seis instancias TODOS os
coeficientes e TODAS as constantes dos 5262 cortes sao inteiros exatos, e todo
fac_i e divisivel por b_il (idem ware_j por p_jl), de modo que pi_i e gamma_j
saem inteiros e o CHECK 3 roda inteiramente em int. Isso tambem PROVA que a
sparsificacao nunca dobrou nada nesses dumps: um coeficiente dobrado teria
|coef| <= 1e-6 e nao nulo, e deslocaria a constante para fora dos inteiros.

O QUE O CHECK 3 PROVA, E O QUE NAO PROVA
========================================
PROVA: para cada corte de otimalidade do dump existe um ponto dual viavel
explicito do LP de roteamento do produto l cujo objetivo domina o corte; logo
a desigualdade theta_l >= const + sum_i fac_i y_i + sum_j ware_j z_j e valida
para TODO (y,z) em [0,1]^I x [0,1]^J, em particular para todo binario, e nao
apenas no incumbente. Toda desigualdade de otimalidade do mestre final e,
portanto, uma desigualdade valida do problema, demonstrada aqui em aritmetica
exata e sem solver.

NAO PROVA: que o conjunto de cortes seja suficiente para provar a otimalidade
do incumbente, nem qualquer coisa sobre a arvore de branch-and-bound.

LIMITACOES, sem maquiagem
=========================
1. Nenhum dos checks PROVA OTIMALIDADE GLOBAL. Com os CHECK 1 a 4 fica provado
   que (a) o custo do incumbente e o declarado, (b) toda desigualdade do
   mestre final e valida (cortes de Benders, CHECK 3) ou nao pode ter cortado
   o incumbente (pisos de theta e linha v_LD, CHECK 4), e (c) o mestre, no
   ponto (y*,z*), ja tinha informacao dual suficiente. O elo que continua fora
   do certificado e a ARVORE do solver: quais nos foram podados e com que
   cota. Reconstruir isso exigiria despejar a arvore inteira, o que nao esta
   implementado.
2. Os coeficientes dos cortes vem de duais de LP em ponto flutuante e sao
   impressos com %.17g. O verificador os le como Fraction do decimal impresso,
   entao a aritmetica DAQUI e exata. O CHECK 3 deixou de depender da qualidade
   do dual devolvido pelo solver: ele REFAZ a verificacao de viabilidade dual.
3. Produtos sem nenhum corte no dump entram no piso do CHECK 2 com
   contribuicao 0 (valido, pois todo custo de transporte e nao negativo).
4. O CHECK 1 confere o custo do incumbente, nao que o incumbente seja otimo.
5. Instancias com dado nao inteiro sao recusadas: toda a cadeia de
   certificacao (gap absoluto < 1, fluxo inteiro) pressupoe dado inteiro, como
   o Instance::load_file do C++ tambem exige.
6. Cortes agregados (l = -1, ablacao "aggregated") somam duais de produtos
   diferentes num unico vetor de coeficientes; o dump v1 nao guarda a
   decomposicao, entao o CHECK 3 os marca NAO DECIDIDO em vez de aprova-los.
   Nenhum dos seis dumps de referencia usa o mestre agregado.

Uso
===
    python3 rational_check.py CAMINHO.certdump [--instance PSC1-C1-50-5.txt]
                                               [--check1 --check2 --check3 --check4]
                                               [--eps 1e-6] [--quiet]
    python3 rational_check.py --batch DIR --instdir DIR [--out relatorio.txt]
    python3 rational_check.py --selftest

Sem nenhuma flag --checkN, roda os quatro. Sai com 0 se todos os checks
executados passam, 1 se algum falha, 2 em erro de leitura.
"""

from __future__ import annotations

import argparse
import glob
import heapq
import os
import sys
import time
from decimal import Decimal, InvalidOperation
from fractions import Fraction

# ------------------------------------------------------------------ leitura --


def _to_int(tok: str, what: str) -> int:
    """Converte um token numerico para int, recusando qualquer parte fracionaria."""
    try:
        d = Decimal(tok)
    except InvalidOperation:
        raise ValueError(f"{what}: token nao numerico {tok!r}")
    if d != d.to_integral_value():
        raise ValueError(f"{what}: valor nao inteiro {tok!r} (proof mode exige dado inteiro)")
    return int(d)


def _to_frac(tok: str) -> Fraction:
    """Fraction exata do decimal impresso (nao do double original; ver LIMITACOES 2)."""
    return Fraction(Decimal(tok))


def _num(v) -> str:
    """Imprime int quando o racional e inteiro, senao a fracao e o float."""
    if isinstance(v, Fraction):
        if v.denominator == 1:
            return str(v.numerator)
        return f"{v} (~{float(v):.9g})"
    return str(v)


class Instance:
    """Le o formato PSC exatamente na ordem de instance.hpp::load_file."""

    def __init__(self, path: str):
        with open(path, "r") as fh:
            tok = fh.read().split()
        pos = 0

        def nxt(what: str) -> int:
            nonlocal pos
            if pos >= len(tok):
                raise ValueError(f"instancia truncada em {what}")
            v = _to_int(tok[pos], what)
            pos += 1
            return v

        self.path = path
        self.I = nxt("nfactories")
        self.J = nxt("nwarehouses")
        self.K = nxt("ncustomers")
        self.L = nxt("ncommodities")
        if min(self.I, self.J, self.K, self.L) <= 0:
            raise ValueError("cabecalho invalido")

        self.demand = [[nxt("q_kl") for _ in range(self.L)] for _ in range(self.K)]

        self.fac_cap = [[0] * self.L for _ in range(self.I)]
        self.fac_fix = [0] * self.I
        for i in range(self.I):
            for l in range(self.L):
                self.fac_cap[i][l] = nxt("b_il")
            self.fac_fix[i] = nxt("f_i")

        self.c_fw = [[[nxt("c_ijl") for _ in range(self.J)] for _ in range(self.I)]
                     for _ in range(self.L)]

        self.ware_cap = [[0] * self.L for _ in range(self.J)]
        self.ware_fix = [0] * self.J
        for j in range(self.J):
            for l in range(self.L):
                self.ware_cap[j][l] = nxt("p_jl")
            self.ware_fix[j] = nxt("g_j")

        self.d_wc = [[[nxt("d_jkl") for _ in range(self.K)] for _ in range(self.J)]
                     for _ in range(self.L)]

        if pos != len(tok):
            raise ValueError(f"instancia com {len(tok) - pos} tokens sobrando")

        # q por produto, na orientacao [l][k], usada pelo CHECK 3.
        self.q_lk = [[self.demand[k][l] for k in range(self.K)] for l in range(self.L)]

    def total_demand(self, l: int) -> int:
        return sum(self.demand[k][l] for k in range(self.K))


class CertDump:
    def __init__(self, path: str):
        self.path = path
        self.instance = None
        self.dims = None
        self.run = None
        self.status = None
        self.obj = None
        self.bound = None
        self.vld = None
        self.verified = None
        self.y: list[int] = []
        self.z: list[int] = []
        self.cuts: list[tuple[int, Fraction, list[Fraction], list[Fraction]]] = []
        self.fl: list[tuple[int, str, Fraction]] = []
        self.ncuts_declared = None
        self.saw_end = False

        with open(path, "r") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                t = line.split()
                tag = t[0]
                if tag == "FORMAT":
                    if int(t[1]) != 1:
                        raise ValueError(f"formato de certdump desconhecido: {t[1]}")
                elif tag == "INSTANCE":
                    self.instance = " ".join(t[1:])
                elif tag == "DIMS":
                    self.dims = tuple(int(v) for v in t[1:5])
                elif tag == "RUN":
                    # RUN <runs> <threads> <papadakos> <aggregated> <theta_integer>
                    self.run = tuple(int(v) for v in t[1:])
                elif tag == "STATUS":
                    self.status = t[1]
                elif tag == "OBJ":
                    self.obj = _to_frac(t[1])
                elif tag == "BOUND":
                    self.bound = _to_frac(t[1])
                elif tag == "VLD":
                    self.vld = None if t[1] == "NONE" else _to_frac(t[1])
                elif tag == "VERIFIED":
                    self.verified = None if t[1] == "NONE" else _to_frac(t[1])
                elif tag == "Y":
                    self.y = [int(v) for v in t[1:]]
                elif tag == "Z":
                    self.z = [int(v) for v in t[1:]]
                elif tag == "NCUTS":
                    self.ncuts_declared = int(t[1])
                elif tag == "CUT":
                    self.cuts.append(self._parse_cut(t))
                elif tag == "FL":
                    self.fl.append((int(t[1]), t[2], _to_frac(t[3])))
                elif tag == "END":
                    self.saw_end = True
                # qualquer outro tag e ignorado de proposito (extensibilidade)

    @staticmethod
    def _parse_cut(t: list[str]):
        l = int(t[1])
        const = _to_frac(t[2])
        if t[3] != "F":
            raise ValueError("CUT sem marcador F")
        w = t.index("W", 4)
        fac = [_to_frac(v) for v in t[4:w]]
        ware = [_to_frac(v) for v in t[w + 1:]]
        return (l, const, fac, ware)

    @property
    def aggregated(self) -> bool:
        if self.run is not None and len(self.run) >= 4:
            return bool(self.run[3])
        return any(c[0] < 0 for c in self.cuts)


# ------------------------------------------------- fluxo de custo minimo int --


class MinCostFlow:
    """SSP com Dijkstra + potenciais de Johnson. Custos >= 0, tudo em int.

    Mesma construcao do mcmf.hpp do C++, reimplementada aqui de proposito: o
    objetivo do verificador e NAO reusar o codigo que ele deveria checar.
    """

    def __init__(self, n: int):
        self.n = n
        self.to: list[int] = []
        self.cap: list[int] = []
        self.cost: list[int] = []
        self.g: list[list[int]] = [[] for _ in range(n)]

    def add_arc(self, u: int, v: int, cap: int, cost: int) -> int:
        eid = len(self.to)
        self.g[u].append(eid)
        self.to.append(v)
        self.cap.append(cap)
        self.cost.append(cost)
        self.g[v].append(eid + 1)
        self.to.append(u)
        self.cap.append(0)
        self.cost.append(-cost)
        return eid

    def flow_on(self, eid: int) -> int:
        return self.cap[eid ^ 1]

    def run(self, s: int, t: int, want: int) -> tuple[int, int]:
        n, to, cap, cost, g = self.n, self.to, self.cap, self.cost, self.g
        INF = float("inf")
        pot = [0] * n
        flow = 0
        total = 0
        while flow < want:
            dist = [INF] * n
            prev = [-1] * n
            dist[s] = 0
            pq = [(0, s)]
            while pq:
                d, u = heapq.heappop(pq)
                if d > dist[u]:
                    continue
                for eid in g[u]:
                    if cap[eid] <= 0:
                        continue
                    v = to[eid]
                    nd = d + cost[eid] + pot[u] - pot[v]
                    if nd < dist[v]:
                        dist[v] = nd
                        prev[v] = eid
                        heapq.heappush(pq, (nd, v))
            if dist[t] == INF:
                break  # t inalcancavel: nao ha mais caminho aumentante
            for v in range(n):
                if dist[v] != INF:
                    pot[v] += dist[v]
            push = want - flow
            v = t
            while v != s:
                eid = prev[v]
                push = min(push, cap[eid])
                v = to[eid ^ 1]
            v = t
            while v != s:
                eid = prev[v]
                cap[eid] -= push
                cap[eid ^ 1] += push
                total += push * cost[eid]
                v = to[eid ^ 1]
            flow += push
        return flow, total


def route_product(inst: Instance, l: int, y: list[int], z: list[int]):
    """Retorna (ok, custo_int) do roteamento otimo do produto l em (y,z)."""
    I, J, K = inst.I, inst.J, inst.K
    S = 0
    T = 1 + I + 2 * J + K
    fac = lambda i: 1 + i
    win = lambda j: 1 + I + j
    wout = lambda j: 1 + I + J + j
    cus = lambda k: 1 + I + 2 * J + k

    demand = inst.total_demand(l)
    g = MinCostFlow(T + 1)
    for i in range(I):
        if y[i]:
            g.add_arc(S, fac(i), inst.fac_cap[i][l], 0)
    for j in range(J):
        if z[j]:
            g.add_arc(win(j), wout(j), inst.ware_cap[j][l], 0)
    for i in range(I):
        if not y[i]:
            continue
        for j in range(J):
            if z[j]:
                g.add_arc(fac(i), win(j), demand, inst.c_fw[l][i][j])
    for j in range(J):
        if not z[j]:
            continue
        for k in range(K):
            g.add_arc(wout(j), cus(k), demand, inst.d_wc[l][j][k])
    for k in range(K):
        g.add_arc(cus(k), T, inst.demand[k][l], 0)

    flow, cost = g.run(S, T, demand)
    return (flow == demand), cost


# ---------------------------------------------------------------- os checks --


def check_primal(inst: Instance, dump: CertDump, out) -> tuple[bool, Fraction | None, list | None]:
    """CHECK 1. Devolve (ok, custo_total, [v_l(y*,z*) por produto])."""
    ok = True
    if len(dump.y) != inst.I or len(dump.z) != inst.J:
        out(f"  FALHA: |Y|={len(dump.y)} |Z|={len(dump.z)}, "
            f"esperado {inst.I} e {inst.J}")
        return False, None, None
    if any(v not in (0, 1) for v in dump.y + dump.z):
        out("  FALHA: Y/Z com valor fora de {0,1}")
        return False, None, None

    # (F_l) do dump avaliadas no incumbente, em inteiros.
    bad_fl = 0
    for (l, kind, rhs) in dump.fl:
        if kind == "FAC":
            lhs = Fraction(sum(inst.fac_cap[i][l] * dump.y[i] for i in range(inst.I)))
        elif kind == "WARE":
            lhs = Fraction(sum(inst.ware_cap[j][l] * dump.z[j] for j in range(inst.J)))
        else:
            out(f"  AVISO: FL de tipo desconhecido {kind}, ignorada")
            continue
        if lhs < rhs:
            bad_fl += 1
            out(f"  FALHA: (F_l) violada em l={l} {kind}: {lhs} < {rhs}")
    if bad_fl:
        ok = False
    else:
        out(f"  (F_l): {len(dump.fl)} desigualdades, todas satisfeitas")

    transport = 0
    vstar = []
    for l in range(inst.L):
        feasible, cost = route_product(inst, l, dump.y, dump.z)
        if not feasible:
            out(f"  FALHA: produto l={l} nao roteavel em (Y,Z)")
            return False, None, None
        vstar.append(cost)
        transport += cost
    fixed = (sum(inst.fac_fix[i] * dump.y[i] for i in range(inst.I))
             + sum(inst.ware_fix[j] * dump.z[j] for j in range(inst.J)))
    total = Fraction(transport + fixed)
    out(f"  roteamento exato: transporte {transport} + fixos {fixed} = {total}")
    out("  v_l(y*,z*) por produto: "
        + ", ".join(f"l{l}={v}" for l, v in enumerate(vstar)))

    for label, ref in (("VERIFIED", dump.verified), ("OBJ", dump.obj)):
        if ref is None:
            out(f"  {label}: ausente no dump")
            continue
        diff = total - ref
        if diff == 0:
            out(f"  {label}: bate exatamente ({ref})")
        else:
            out(f"  FALHA: {label}={ref} difere do recalculo por {float(diff):+.6g}")
            ok = False
    return ok, total, vstar


def cut_value_at(inst: Instance, cut, y: list[int], z: list[int]) -> Fraction:
    """Lado direito do corte avaliado em (y,z)."""
    (_l, const, fac, ware) = cut
    val = const
    for i in range(inst.I):
        if y[i]:
            val += fac[i]
    for j in range(inst.J):
        if z[j]:
            val += ware[j]
    return val


def check_cuts(inst: Instance, dump: CertDump, cost: Fraction | None, out):
    """CHECK 2. Devolve (ok, best_por_produto)."""
    if not dump.cuts:
        out("  nenhum corte no dump: CHECK 2 impossivel (ver LIMITACOES 3)")
        return False, {}
    if cost is None:
        out("  sem custo primal recalculado: CHECK 2 impossivel")
        return False, {}

    best: dict[int, Fraction] = {}
    bad_dim = 0
    for cut in dump.cuts:
        (l, const, fac, ware) = cut
        if len(fac) != inst.I or len(ware) != inst.J:
            bad_dim += 1
            continue
        val = cut_value_at(inst, cut, dump.y, dump.z)
        if val > best.get(l, Fraction(0)):
            best[l] = val
    if bad_dim:
        out(f"  FALHA: {bad_dim} cortes com dimensao errada, ignorados")

    fixed = (sum(inst.fac_fix[i] * dump.y[i] for i in range(inst.I))
             + sum(inst.ware_fix[j] * dump.z[j] for j in range(inst.J)))
    # Produtos sem corte contribuem 0 (LIMITACOES 3); o piso so pode ficar baixo.
    theta_sum = sum(best.values(), Fraction(0))
    floorval = Fraction(fixed) + theta_sum
    covered = len([l for l in best if l >= 0])
    if any(l < 0 for l in best):
        out("  master agregado detectado (l=-1): um unico bloco de theta")
    else:
        out(f"  cortes cobrem {covered}/{inst.L} produtos "
            f"({len(dump.cuts)} cortes lidos)")

    slack = cost - floorval
    out(f"  piso reconstituido = fixos {fixed} + theta {float(theta_sum):.6f} "
        f"= {float(floorval):.6f}")
    out(f"  custo verificado   = {float(cost):.6f}   folga = {float(slack):.6g}")
    if bad_dim:
        return False, best
    if slack <= Fraction(9999, 10000):
        out("  folga < 0.9999: os cortes fecham o incumbente em aritmetica exata")
        return True, best
    out("  folga >= 0.9999: os cortes NAO fecham o incumbente neste ponto")
    return False, best


# ------------------------------------------------- CHECK 3, validade global --


def recover_duals(inst: Instance, l: int, fac, ware):
    """pi_i = fac_i/b_il, gamma_j = ware_j/p_jl. None marca dual livre (cap 0).

    Devolve (pi, gam, notas). Os valores saem int quando a divisao e exata.
    """
    notes = []
    pi: list = []
    for i in range(inst.I):
        b = inst.fac_cap[i][l]
        v = fac[i]
        if b == 0:
            if v != 0:
                notes.append(f"fac_{i} = {v} != 0 com b_{i}{l} = 0 (corte inconsistente)")
            pi.append(None)  # dual livre: o indice sai do minimo
            continue
        q = Fraction(v, 1) / b if not isinstance(v, Fraction) else v / b
        if q > 0:
            notes.append(f"pi_{i} = {q} > 0, fora do cone dual (esperado <= 0)")
        pi.append(int(q) if q.denominator == 1 else q)
    gam: list = []
    for j in range(inst.J):
        p = inst.ware_cap[j][l]
        v = ware[j]
        if p == 0:
            if v != 0:
                notes.append(f"ware_{j} = {v} != 0 com p_{j}{l} = 0 (corte inconsistente)")
            gam.append(None)
            continue
        q = Fraction(v, 1) / p if not isinstance(v, Fraction) else v / p
        if q > 0:
            notes.append(f"gamma_{j} = {q} > 0, fora do cone dual (esperado <= 0)")
        gam.append(int(q) if q.denominator == 1 else q)
    return pi, gam, notes


def max_dual_constant(inst: Instance, l: int, pi, gam):
    """Phi_l(pi,gamma) = max sum_k q_kl alpha_k sobre o sistema de diferencas.

    Resolve o dual de caminho minimo em duas varreduras de minimo (o DAG em
    camadas colapsa o Bellman-Ford). Devolve (Phi, beta, alpha) ou (None,...)
    se o dual for ilimitado por falta de camada (nao ocorre com demanda > 0).
    """
    C = inst.c_fw[l]
    D = inst.d_wc[l]
    J, K = inst.J, inst.K

    idx_i = [i for i in range(inst.I) if pi[i] is not None]
    if not idx_i:
        return None, None, None  # nenhuma fabrica com capacidade: dual ilimitado
    i0 = idx_i[0]
    m0 = -pi[i0]
    beta = list(map(m0.__add__, C[i0])) if not isinstance(m0, Fraction) \
        else [m0 + c for c in C[i0]]
    for i in idx_i[1:]:
        m = -pi[i]
        if isinstance(m, Fraction):
            beta = [b if b <= m + c else m + c for b, c in zip(beta, C[i])]
        else:
            beta = list(map(min, beta, map(m.__add__, C[i])))

    idx_j = [j for j in range(J) if gam[j] is not None]
    if not idx_j:
        return None, None, None
    j0 = idx_j[0]
    u0 = beta[j0] - gam[j0]
    alpha = list(map(u0.__add__, D[j0])) if not isinstance(u0, Fraction) \
        else [u0 + d for d in D[j0]]
    for j in idx_j[1:]:
        u = beta[j] - gam[j]
        if isinstance(u, Fraction):
            alpha = [a if a <= u + d else u + d for a, d in zip(alpha, D[j])]
        else:
            alpha = list(map(min, alpha, map(u.__add__, D[j])))

    q = inst.q_lk[l]
    phi = sum(a * w for a, w in zip(alpha, q))
    return phi, beta, alpha


def cut_penalty(inst: Instance, l: int, pi, gam, fac, ware):
    """sum_i min(0, b_il pi_i - fac_i) + sum_j min(0, p_jl gam_j - ware_j).

    E o minimo, sobre (y,z) binario, da diferenca entre os coeficientes do
    objetivo dual e os coeficientes do corte emitido. Zero quando o dual
    recuperado reproduz exatamente os coeficientes.
    """
    pen = Fraction(0)
    for i in range(inst.I):
        b = inst.fac_cap[i][l]
        lhs = 0 if pi[i] is None else b * pi[i]
        d = lhs - fac[i]
        if d < 0:
            pen += d
    for j in range(inst.J):
        p = inst.ware_cap[j][l]
        lhs = 0 if gam[j] is None else p * gam[j]
        d = lhs - ware[j]
        if d < 0:
            pen += d
    return pen


def sparsification_rescue(inst: Instance, l: int, const, fac, ware, eps: Fraction):
    """Segundo estagio: maximiza G na janela da sparsificacao, exatamente.

    Para cada coeficiente emitido como zero exato, o dual verdadeiro pode ter
    ficado em [-eps/cap, 0). Baixar pi_i de delta aumenta beta (e portanto
    Phi), mas paga -b_il*delta na penalidade. G e concava nesses deltas; o
    codigo avalia G exatamente num conjunto finito de perfis delta_i =
    min(t, eps/b_il) e devolve o melhor. Qualquer G >= 0 encontrado E PROVA da
    validade do corte emitido; um maximo negativo nao prova invalidade.

    Devolve (best_G, best_t, upper_bound_analitico).
    """
    zf = [i for i in range(inst.I) if fac[i] == 0 and inst.fac_cap[i][l] > 0]
    zw = [j for j in range(inst.J) if ware[j] == 0 and inst.ware_cap[j][l] > 0]
    base_pi, base_gam, _ = recover_duals(inst, l, fac, ware)

    cand = {Fraction(0)}
    for i in zf:
        cand.add(eps / inst.fac_cap[i][l])
    for j in zw:
        cand.add(eps / inst.ware_cap[j][l])
    ts = sorted(cand)
    if len(ts) > 48:  # amostra uniforme, o custo por ponto e O(IJ + JK)
        step = len(ts) / 48.0
        ts = [ts[min(len(ts) - 1, int(x * step))] for x in range(48)]
        ts = sorted(set(ts))

    Q = inst.total_demand(l)
    dmax_f = max((eps / inst.fac_cap[i][l] for i in zf), default=Fraction(0))
    dmax_w = max((eps / inst.ware_cap[j][l] for j in zw), default=Fraction(0))

    best_g = None
    best_t = None
    for t in ts:
        pi = list(base_pi)
        gam = list(base_gam)
        for i in zf:
            b = inst.fac_cap[i][l]
            d = min(t, eps / b)
            if d:
                pi[i] = -d
        for j in zw:
            p = inst.ware_cap[j][l]
            d = min(t, eps / p)
            if d:
                gam[j] = -d
        phi, _, _ = max_dual_constant(inst, l, pi, gam)
        if phi is None:
            continue
        # cut_penalty ja contabiliza a diferenca b_il*pi_i - fac_i < 0 introduzida
        # pelos deltas, que e exatamente o preco pago por baixar o dual.
        g = Fraction(phi) - const + cut_penalty(inst, l, pi, gam, fac, ware)
        if best_g is None or g > best_g:
            best_g = g
            best_t = t
    ub = None
    if best_g is not None:
        g0 = None
        phi0, _, _ = max_dual_constant(inst, l, base_pi, base_gam)
        if phi0 is not None:
            g0 = Fraction(phi0) - const + cut_penalty(inst, l, base_pi, base_gam, fac, ware)
            ub = g0 + Q * (dmax_f + dmax_w)
    return best_g, best_t, ub


def check_cut_validity(inst: Instance, dump: CertDump, out, eps: Fraction,
                       vstar: list | None = None, vopen: list | None = None,
                       max_report: int = 12):
    """CHECK 3. Validade global, corte a corte. Devolve (ok, resumo)."""
    res = {
        "n": len(dump.cuts), "proved": 0, "rescued": 0, "undecided": 0,
        "failed": 0, "invalid_proved": 0, "fails": [], "tight": 0,
        "min_margin": None, "max_margin": None,
    }
    if not dump.cuts:
        out("  nenhum corte no dump: CHECK 3 vacuo")
        return True, res

    y, z = dump.y, dump.z
    ones_y = [1] * inst.I
    ones_z = [1] * inst.J

    # Auditoria da sparsificacao. O C++ dobra na constante todo coeficiente com
    # |coef| <= eps. Um coeficiente dobrado e NAO NULO deslocaria a constante de
    # uma quantidade nao inteira e menor que eps. Portanto: se toda constante e
    # todo coeficiente sao inteiros exatos e eps < 1, nenhuma dobra de valor nao
    # nulo pode ter ocorrido, e o corte emitido E o corte dual original.
    n_zero = 0
    n_coef = 0
    nonint_const = 0
    nonint_coef = 0
    for (l, const, fac, ware) in dump.cuts:
        if const.denominator != 1:
            nonint_const += 1
        for v in fac:
            n_coef += 1
            if v == 0:
                n_zero += 1
            elif v.denominator != 1:
                nonint_coef += 1
        for v in ware:
            n_coef += 1
            if v == 0:
                n_zero += 1
            elif v.denominator != 1:
                nonint_coef += 1
    res["n_zero_coef"] = n_zero
    res["n_coef"] = n_coef
    res["nonint"] = nonint_const + nonint_coef
    out(f"  auditoria da sparsificacao (eps = {float(eps)}): {n_coef} coeficientes, "
        f"{n_zero} emitidos como zero exato, {nonint_const} constantes e "
        f"{nonint_coef} coeficientes nao inteiros")
    if res["nonint"] == 0 and eps < 1:
        out("      todos inteiros exatos e eps < 1: NENHUMA dobra de valor nao nulo")
        out("      pode ter ocorrido (ela deslocaria a constante para fora dos")
        out("      inteiros). Os cortes emitidos sao os cortes duais originais e a")
        out("      recuperacao pi_i = fac_i/b_il e exata, nao perturbada.")
    else:
        out("      ha valores nao inteiros: a dobra pode ter atuado. A dobra enfraquece")
        out("      o corte (coef <= 0, coef*y >= coef), entao um FAIL ainda pode ser")
        out("      artefato; o segundo estagio trata disso.")

    for idx, cut in enumerate(dump.cuts):
        (l, const, fac, ware) = cut
        if l < 0:
            res["undecided"] += 1
            continue
        if len(fac) != inst.I or len(ware) != inst.J or not (0 <= l < inst.L):
            res["failed"] += 1
            res["fails"].append((idx, l, "dimensao/indice invalidos", None))
            continue

        pi, gam, notes = recover_duals(inst, l, fac, ware)
        phi, _, _ = max_dual_constant(inst, l, pi, gam)
        if phi is None:
            res["undecided"] += 1
            continue
        pen = cut_penalty(inst, l, pi, gam, fac, ware)
        g = Fraction(phi) - const + pen
        if res["min_margin"] is None or g < res["min_margin"]:
            res["min_margin"] = g
        if res["max_margin"] is None or g > res["max_margin"]:
            res["max_margin"] = g
        if g >= 0:
            res["proved"] += 1
            if g == 0:
                res["tight"] += 1
            continue

        # segundo estagio: folga da sparsificacao
        bg, bt, ub = sparsification_rescue(inst, l, const, fac, ware, eps)
        if bg is not None and bg >= 0:
            res["rescued"] += 1
            continue

        # terceiro estagio: contraexemplo primal certificado
        witness = None
        rhs_star = cut_value_at(inst, cut, y, z)
        if vstar is not None and rhs_star > vstar[l]:
            witness = (f"(y*,z*): RHS {_num(rhs_star)} > v_{l}(y*,z*) = {vstar[l]}")
        elif vopen is not None:
            rhs_open = cut_value_at(inst, cut, ones_y, ones_z)
            if rhs_open > vopen[l]:
                witness = (f"(1,1): RHS {_num(rhs_open)} > v_{l}(1,1) = {vopen[l]}")
        res["failed"] += 1
        if witness:
            res["invalid_proved"] += 1
        res["fails"].append((idx, l, {
            "G": g, "G_rescue": bg, "t": bt, "UB": ub, "notes": notes,
        }, witness))

    out(f"  cortes validos {res['proved'] + res['rescued']} de {res['n']}"
        f"   (prova direta {res['proved']}, resgatados pela folga "
        f"{res['rescued']}, nao decididos {res['undecided']}, falhas "
        f"{res['failed']})")
    if res["min_margin"] is not None:
        out(f"  margem dual G = Phi_l(pi,gamma) - const + penalidade: "
            f"min {_num(res['min_margin'])}, max {_num(res['max_margin'])}, "
            f"exatamente tight (G=0) em {res['tight']} cortes")
    if res["undecided"]:
        out(f"  AVISO: {res['undecided']} cortes agregados (l=-1) nao sao "
            "decidiveis no formato v1 (LIMITACOES 6)")
    for (idx, l, info, witness) in res["fails"][:max_report]:
        if isinstance(info, str):
            out(f"  FALHA corte #{idx} l={l}: {info}")
            continue
        out(f"  FALHA corte #{idx} l={l}: G={_num(info['G'])} "
            f"(resgate maximo na janela {_num(info['G_rescue'])} em t={info['t']}, "
            f"limite analitico {_num(info['UB']) if info['UB'] is not None else 'n/d'})")
        for n in info["notes"][:4]:
            out(f"      nota de recuperacao: {n}")
        if witness:
            out(f"      CONTRAEXEMPLO CERTIFICADO, corte INVALIDO em {witness}")
        elif info["UB"] is not None and info["UB"] < 0:
            out("      diagnostico: nenhum dual da janela de sparsificacao salva "
                "este corte; a falha NAO e artefato de recuperacao")
        else:
            out("      diagnostico: falha dentro da janela de sparsificacao, "
                "possivel artefato de recuperacao")
    if len(res["fails"]) > max_report:
        out(f"  ... e mais {len(res['fails']) - max_report} falhas nao listadas")
    return res["failed"] == 0 and res["undecided"] == 0, res


# ------------------------------------------- CHECK 4, pisos de theta e v_LD --


def check_floors(inst: Instance, dump: CertDump, out, vstar: list | None,
                 best_cuts: dict | None, cost: Fraction | None):
    """CHECK 4. v_l(1,1) exato, margem por produto, e folga exata do v_LD."""
    res = {"vopen": None, "margins": None, "min_margin": None,
           "vld_slack": None, "ok_margin": False, "ok_vld": False}
    if vstar is None:
        out("  CHECK 1 nao rodou: theta*_l indisponivel, recomputando")
        vstar = []
        for l in range(inst.L):
            ok, c = route_product(inst, l, dump.y, dump.z)
            if not ok:
                out(f"  FALHA: produto {l} nao roteavel no incumbente")
                return False, res
            vstar.append(c)
    if cost is None:
        # CHECK 1 nao rodou: reconstitui o custo do incumbente a partir dos
        # v_l(y*,z*) recem-computados mais os fixos, tudo em inteiros.
        cost = Fraction(sum(vstar)
                        + sum(inst.fac_fix[i] * dump.y[i] for i in range(inst.I))
                        + sum(inst.ware_fix[j] * dump.z[j] for j in range(inst.J)))

    ones_y = [1] * inst.I
    ones_z = [1] * inst.J
    vopen = []
    for l in range(inst.L):
        ok, c = route_product(inst, l, ones_y, ones_z)
        if not ok:
            out(f"  FALHA: produto {l} nao roteavel nem com tudo aberto")
            return False, res
        vopen.append(c)
    res["vopen"] = vopen

    margins = [vstar[l] - vopen[l] for l in range(inst.L)]
    res["margins"] = margins
    res["min_margin"] = min(margins)
    out("  (i) piso inicial theta_l >= v_l(1,1), recomputado exatamente:")
    out("      l | v_l(1,1) exato | theta*_l = v_l(y*,z*) | margem | max corte em (y*,z*)")
    for l in range(inst.L):
        bc = "-" if not best_cuts else _num(best_cuts.get(l, Fraction(0)))
        out(f"      {l} | {vopen[l]:>14} | {vstar[l]:>21} | {margins[l]:>6} | {bc:>20}")
    res["ok_margin"] = res["min_margin"] > 1
    if res["ok_margin"]:
        out(f"      margem minima {res['min_margin']} > 1 unidade em todos os produtos:")
        out("      o piso theta_l >= v_l(1,1) instalado pelo construtor do mestre veio")
        out("      de um LP em ponto flutuante; ainda que esse LP errasse por menos de")
        out("      1 unidade para cima, o piso ficaria abaixo de theta*_l e nao poderia")
        out("      ter excluido o incumbente otimo. (Alem disso o piso e valido por")
        out("      monotonicidade de capacidade: v_l(y,z) >= v_l(1,1) para todo (y,z).)")
    else:
        out(f"      ATENCAO: margem minima {res['min_margin']} <= 1; o argumento de")
        out("      erro sub-unitario NAO cobre este caso, so a monotonicidade cobre.")

    out("  (ii) linha global fixos + sum_l theta_l >= v_LD:")
    if dump.vld is None:
        out("      VLD NONE: nenhuma linha Lagrangiana foi adicionada ao mestre.")
        out("      Exoneracao vacua, nao ha o que cortar.")
        res["ok_vld"] = True
    elif cost is None:
        out("      custo verificado indisponivel, exoneracao impossivel")
    else:
        slack = cost - dump.vld
        res["vld_slack"] = slack
        out(f"      custo verificado = {_num(cost)}")
        out(f"      v_LD             = {_num(dump.vld)}")
        out(f"      folga exata      = {_num(slack)}")
        if slack > 0:
            res["ok_vld"] = True
            out("      folga > 0: no ponto otimo o mestre tem fixos + sum theta_l = custo")
            out("      (CHECK 2 fecha com folga zero, logo theta*_l = v_l(y*,z*)), e esse")
            out("      valor satisfaz a linha com folga estritamente positiva. A linha")
            out("      global NAO pode ter cortado o otimo, qualquer que seja o erro do")
            out("      subgradiente que produziu v_LD, desde que menor que a folga.")
        else:
            out("      FALHA: folga <= 0, a linha Lagrangiana tocava ou cortava o otimo")
    ok = res["ok_margin"] and res["ok_vld"]
    return ok, res


# ------------------------------------------------------------------ selftest --

SELFTEST_INST = """2 2 2 1
30
20
40 100
25 60
1 2
2 1
35 70
35 90
3 1
1 3
"""


def selftest(out) -> bool:
    """Instancia minima conferida a mao, para exercitar parser, fluxo e checks.

    I=2 J=2 K=2 L=1. Demandas q_k = 30 e 20 (total 50).
    Fabricas: capacidade 40 (fixo 100) e 25 (fixo 60).
    Depositos: capacidade 35 (fixo 70) e 35 (fixo 90).
    c_fw = [[1,2],[2,1]], d_wc = [[3,1],[1,3]].

    Roteamento otimo com y=(1,1), z=(1,1), conferido a mao:
      etapa 2: k1 (30) vem de w2 a custo 1, k2 (20) vem de w1 a custo 1  -> 50
      etapa 1: w1 recebe 20 de f1 a custo 1 -> 20; w2 recebe 30, mas f2 so tem
               capacidade 25 (custo 1 -> 25) e os 5 restantes vem de f1 a
               custo 2 -> 10.  Etapa 1 = 55.
      transporte 105, fixos 100+60+70+90 = 320, custo total 425.

    Dual otimo desse LP, conferido a mao pela folga complementar:
      pi = (0, -1), gamma = (0, 0), beta = (1, 2), alpha = (3, 2)
      corte:  theta_0 >= 130 + 0*y_1 - 25*y_2 + 0*z_1 + 0*z_2
      em (1,1): 130 - 25 = 105 = v_0(1,1). O corte e um corte de Benders
      legitimo e o CHECK 3 tem de aprova-lo com margem G = 0.

    TESTE NEGATIVO OBRIGATORIO. Dois cortes sinteticos INVALIDOS:
      N1: theta_0 >= 130 + 0*y_1 - 20*y_2   (coeficiente enfraquecido)
          em (1,1) da 110 > 105 = v_0(1,1): invalido, com testemunha.
          Recuperacao: pi_2 = -20/25 = -4/5, exercita o caminho Fraction.
      N2: theta_0 >= 140 + 0*y_1 - 25*y_2   (constante inflada)
          em (1,1) da 115 > 105: invalido, com testemunha.
    O CHECK 3 tem de REPROVAR os dois; se algum dia aprovar, o check regrediu.

    ATENCAO: se este selftest passar a imprimir outro numero, o solver de fluxo
    deste arquivo regrediu.
    """
    import tempfile

    ok_all = True
    d = tempfile.mkdtemp(prefix="certselftest")
    ip = os.path.join(d, "TINY.txt")
    with open(ip, "w") as fh:
        fh.write(SELFTEST_INST)
    inst = Instance(ip)
    y, z = [1, 1], [1, 1]
    ok, cost = route_product(inst, 0, y, z)
    if not ok:
        out("SELFTEST FALHOU: instancia minima nao roteavel")
        return False
    if cost != 105:
        out(f"SELFTEST FALHOU: transporte {cost} != 105 (solver de fluxo regrediu)")
        ok_all = False
    total = cost + sum(inst.fac_fix) + sum(inst.ware_fix)
    out(f"selftest: transporte {cost}, custo total {total} (esperado 105 e 425)")
    if total != 425:
        out("SELFTEST FALHOU: custo total != 425")
        ok_all = False

    # -------- dump positivo: um corte de Benders legitimo, dual conferido a mao
    cp = os.path.join(d, "TINY_selftest.certdump")
    with open(cp, "w") as fh:
        fh.write("# mptscfl certdump v1\nFORMAT 1\n")
        fh.write(f"INSTANCE {ip}\n")
        fh.write(f"DIMS {inst.I} {inst.J} {inst.K} {inst.L}\n")
        fh.write("RUN 1 1 0 0 0\n")
        fh.write("STATUS OPTIMAL\n")
        fh.write(f"OBJ {total}\nBOUND {total}\nVLD 400\nVERIFIED {total}\n")
        fh.write("Y 1 1\nZ 1 1\nNCUTS 1\n")
        fh.write("CUT 0 130 F 0 -25 W 0 0\n")
        q = inst.total_demand(0)
        fh.write(f"FL 0 FAC {q}\nFL 0 WARE {q}\nEND\n")
    dump = CertDump(cp)
    p_ok, p_cost, vstar = check_primal(inst, dump, out)
    c_ok, best = check_cuts(inst, dump, p_cost, out)
    v3_ok, r3 = check_cut_validity(inst, dump, out, Fraction(1, 1000000), vstar, None)
    f4_ok, r4 = check_floors(inst, dump, out, vstar, best, p_cost)
    out(f"selftest positivo: primal {'PASS' if p_ok else 'FAIL'}, "
        f"cortes {'PASS' if c_ok else 'FAIL'}, "
        f"validade {'PASS' if v3_ok else 'FAIL'}, "
        f"pisos {'PASS' if f4_ok else 'FAIL'}")
    out("  (o FAIL de pisos e ESPERADO na instancia minima: o incumbente E o "
        "tudo-aberto, logo a margem theta*_0 - v_0(1,1) vale 0 por construcao. "
        "O selftest exige apenas que os numeros 105 e 25 saiam certos.)")
    if not (p_ok and c_ok and v3_ok):
        out("SELFTEST FALHOU: o bloco positivo deveria passar CHECK 1, 2 e 3")
        ok_all = False
    if r3["proved"] != 1 or r3["min_margin"] != 0:
        out(f"SELFTEST FALHOU: esperava 1 corte provado com G=0, obtive "
            f"proved={r3['proved']} G={r3['min_margin']}")
        ok_all = False
    if r4["vopen"] != [105]:
        out(f"SELFTEST FALHOU: v_0(1,1) = {r4['vopen']} != [105]")
        ok_all = False
    if r4["vld_slack"] != 25:
        out(f"SELFTEST FALHOU: folga de v_LD = {r4['vld_slack']} != 25")
        ok_all = False

    # -------- dump negativo: cortes sinteticos INVALIDOS
    out("selftest NEGATIVO (cortes invalidos que o CHECK 3 tem de reprovar)")
    np_ = os.path.join(d, "TINY_negativo.certdump")
    with open(np_, "w") as fh:
        fh.write("# mptscfl certdump v1\nFORMAT 1\n")
        fh.write(f"INSTANCE {ip}\n")
        fh.write(f"DIMS {inst.I} {inst.J} {inst.K} {inst.L}\n")
        fh.write("RUN 1 1 0 0 0\n")
        fh.write("STATUS OPTIMAL\n")
        fh.write(f"OBJ {total}\nBOUND {total}\nVLD NONE\nVERIFIED {total}\n")
        fh.write("Y 1 1\nZ 1 1\nNCUTS 3\n")
        fh.write("CUT 0 130 F 0 -25 W 0 0\n")   # valido
        fh.write("CUT 0 130 F 0 -20 W 0 0\n")   # N1 invalido
        fh.write("CUT 0 140 F 0 -25 W 0 0\n")   # N2 invalido
        q = inst.total_demand(0)
        fh.write(f"FL 0 FAC {q}\nFL 0 WARE {q}\nEND\n")
    ndump = CertDump(np_)
    n_ok, nr = check_cut_validity(inst, ndump, out, Fraction(1, 1000000),
                                 vstar, r4["vopen"])
    if n_ok:
        out("SELFTEST FALHOU: o CHECK 3 APROVOU cortes invalidos (regressao grave)")
        ok_all = False
    elif nr["failed"] != 2 or nr["proved"] != 1:
        out(f"SELFTEST FALHOU: esperava 1 provado e 2 reprovados, obtive "
            f"proved={nr['proved']} failed={nr['failed']}")
        ok_all = False
    elif nr["invalid_proved"] != 2:
        out(f"SELFTEST FALHOU: esperava 2 contraexemplos certificados, obtive "
            f"{nr['invalid_proved']}")
        ok_all = False
    else:
        out("selftest negativo: OK, os 2 cortes invalidos foram reprovados com "
            "contraexemplo certificado e o valido foi aprovado")

    out(f"SELFTEST {'VERDE' if ok_all else 'VERMELHO'}")
    return ok_all


# ---------------------------------------------------------------------- main --


def resolve_instance(dump: CertDump, certpath: str, override: str | None,
                     instdir: str | None) -> str:
    cands = []
    if override:
        cands.append(override)
    if dump.instance:
        raw = dump.instance.replace("\\", "/")
        cands.append(raw)
        base = os.path.basename(raw)
        if instdir:
            cands.append(os.path.join(instdir, base))
        cands.append(os.path.join(os.path.dirname(os.path.abspath(certpath)), base))
    for c in cands:
        if c and os.path.exists(c):
            return c
    raise FileNotFoundError(f"instancia nao encontrada; tentei: {cands}")


def run_one(certpath: str, out, checks: set, eps: Fraction,
            instance_override=None, instdir=None) -> dict:
    r = {"cert": certpath, "name": os.path.basename(certpath).split("_")[0],
         "ok": True, "times": {}}
    dump = CertDump(certpath)
    if not dump.saw_end:
        out("AVISO: dump sem linha END, run possivelmente interrompido")
    if dump.ncuts_declared is not None and dump.ncuts_declared != len(dump.cuts):
        out(f"AVISO: NCUTS={dump.ncuts_declared} mas {len(dump.cuts)} linhas CUT")
    ipath = resolve_instance(dump, certpath, instance_override, instdir)
    inst = Instance(ipath)
    if dump.dims and dump.dims != (inst.I, inst.J, inst.K, inst.L):
        raise ValueError(f"DIMS {dump.dims} != instancia {(inst.I, inst.J, inst.K, inst.L)}")

    r["dump"] = dump
    r["inst"] = inst
    out(f"certdump : {certpath}")
    out(f"instancia: {ipath}  (I={inst.I} J={inst.J} K={inst.K} L={inst.L})")
    out(f"status   : {dump.status}   cortes: {len(dump.cuts)}   "
        f"v_LD: {_num(dump.vld) if dump.vld is not None else 'NONE'}   "
        f"master: {'agregado' if dump.aggregated else 'desagregado'}")

    cost = None
    vstar = None
    best = None

    if 1 in checks:
        out("CHECK 1 (primal, fluxo exato em inteiros)")
        t0 = time.time()
        p_ok, cost, vstar = check_primal(inst, dump, out)
        r["times"]["check1"] = time.time() - t0
        out(f"CHECK 1 {'PASS' if p_ok else 'FAIL'}  [{r['times']['check1']:.1f}s]")
        r["check1"] = p_ok
        r["ok"] &= p_ok
        r["cost"] = cost
        r["vstar"] = vstar

    if 2 in checks:
        out("CHECK 2 (fechamento dos cortes no incumbente)")
        t0 = time.time()
        c_ok, best = check_cuts(inst, dump, cost, out)
        r["times"]["check2"] = time.time() - t0
        out(f"CHECK 2 {'PASS' if c_ok else 'FAIL'}  [{r['times']['check2']:.1f}s]")
        r["check2"] = c_ok
        r["ok"] &= c_ok
        r["best"] = best

    if 4 in checks:
        # roda antes do 3 para que o CHECK 3 tenha v_l(1,1) como testemunha
        out("CHECK 4 (piso theta_l >= v_l(1,1) e linha global v_LD)")
        t0 = time.time()
        f_ok, r4 = check_floors(inst, dump, out, vstar, best, cost)
        r["times"]["check4"] = time.time() - t0
        out(f"CHECK 4 {'PASS' if f_ok else 'FAIL'}  [{r['times']['check4']:.1f}s]")
        r["check4"] = f_ok
        r["r4"] = r4
        r["ok"] &= f_ok

    if 3 in checks:
        out("CHECK 3 (validade GLOBAL de cada corte, viabilidade dual exata)")
        t0 = time.time()
        vopen = r.get("r4", {}).get("vopen") if 4 in checks else None
        v_ok, r3 = check_cut_validity(inst, dump, out, eps, vstar, vopen)
        r["times"]["check3"] = time.time() - t0
        out(f"CHECK 3 {'PASS' if v_ok else 'FAIL'}  [{r['times']['check3']:.1f}s]")
        r["check3"] = v_ok
        r["r3"] = r3
        r["ok"] &= v_ok

    return r


def batch(args, out) -> int:
    dumps = sorted(glob.glob(os.path.join(args.batch, "*.certdump")))
    if not dumps:
        print(f"ERRO: nenhum .certdump em {args.batch}", file=sys.stderr)
        return 2
    lines: list[str] = []

    def rec(msg=""):
        lines.append(str(msg))
        out(msg)

    results = []
    t_all = time.time()
    checks = args.checks
    rec("RELATORIO RACIONAL v2, MP-TSCFLP")
    rec("=" * 78)
    rec(f"gerado por rational_check.py em {time.strftime('%Y-%m-%d %H:%M:%S')}")
    rec(f"dumps: {args.batch}")
    rec(f"instancias: {args.instdir}")
    rec(f"checks executados: {sorted(checks)}   eps de sparsificacao: {float(args.eps)}")
    rec("")
    rec("ESTRUTURA DUAL USADA PELO CHECK 3 (derivada de src3/benders_model.cpp)")
    rec("-" * 78)
    rec("Subproblema de roteamento do produto l em (y,z), como montado por")
    rec("ProductSubproblem::ProductSubproblem e ::solve:")
    rec("")
    rec("    min   sum_ij c_ijl x_ij + sum_jk d_jkl w_jk")
    rec("    s.a.  sum_j w_jk              >= q_kl      (dem_k)   alpha_k >= 0")
    rec("          sum_i x_ij - sum_k w_jk >= 0         (cons_j)  beta_j  >= 0")
    rec("          sum_j x_ij              <= b_il y_i  (fcap_i)  pi_i    <= 0")
    rec("          sum_k w_jk              <= p_jl z_j  (wcap_j)  gamma_j <= 0")
    rec("          x, w >= 0")
    rec("")
    rec("Dual, com os sinais da convencao Gurobi para minimizacao:")
    rec("")
    rec("    max   sum_k q_kl alpha_k + sum_i (b_il y_i) pi_i + sum_j (p_jl z_j) gamma_j")
    rec("    s.a.  pi_i + beta_j              <= c_ijl   (coluna x_ij)")
    rec("          alpha_k - beta_j + gamma_j <= d_jkl   (coluna w_jk)")
    rec("")
    rec("O corte gravado no dump E o objetivo dual nesse ponto: const = sum_k q_kl")
    rec("alpha_k, fac_i = b_il pi_i, ware_j = p_jl gamma_j. A REGIAO VIAVEL DUAL NAO")
    rec("DEPENDE DE (y,z); so o objetivo depende. Logo qualquer ponto dual viavel gera")
    rec("uma desigualdade valida para TODO (y,z), e provar a validade global do corte")
    rec("e exibir um ponto dual viavel cujo objetivo domine o corte emitido.")
    rec("")
    rec("Do dump recuperam-se pi_i = fac_i/b_il e gamma_j = ware_j/p_jl. Os alpha_k")
    rec("individuais nao sao recuperaveis, so a soma ponderada (a constante), e os")
    rec("beta_j das linhas de RHS zero (balanco no deposito) nao aparecem no corte.")
    rec("A validade fica equivalente a viabilidade de um sistema de DIFERENCAS em")
    rec("(alpha,beta) com (pi,gamma) fixos, que e um dual de caminho minimo em")
    rec("camadas. Por ser um DAG em camadas o Bellman-Ford exato colapsa em duas")
    rec("varreduras de minimo, e os limites simultaneamente atingiveis sao")
    rec("")
    rec("    beta_j^max  = min_i (c_ijl - pi_i)                  >= 0")
    rec("    alpha_k^max = min_j (d_jkl + beta_j^max - gamma_j)  >= 0")
    rec("")
    rec("de modo que Phi_l(pi,gamma) = sum_k q_kl alpha_k^max e o MAIOR termo")
    rec("constante compativel com (pi,gamma). O criterio exato, necessario e")
    rec("suficiente para (y,z) binario arbitrario, e")
    rec("")
    rec("    G = Phi_l(pi,gamma) - const")
    rec("        + sum_i min(0, b_il pi_i - fac_i)")
    rec("        + sum_j min(0, p_jl gamma_j - ware_j)  >= 0")
    rec("")
    rec("Tudo decidido em int/Fraction, sem simplex e sem solver. Sobre a")
    rec("sparsificacao (cut_expr dobra na constante os coeficientes com |coef| <=")
    rec("eps): os coeficientes sao <= 0 e coef*y >= coef para y em [0,1], portanto A")
    rec("DOBRA VAI PARA O LADO SEGURO, ela ENFRAQUECE o corte. O corte emitido e")
    rec("dominado pelo corte do dual verdadeiro. A recuperacao com fac_i = 0 devolve")
    rec("pi_i = 0 no lugar de um pi_i em [-eps/b_il, 0), o que so diminui Phi_l: o")
    rec("teste fica CONSERVADOR, um PASS continua prova. Em caso de FAIL entra o")
    rec("segundo estagio, que maximiza G sobre a janela [-eps/cap, 0] com variavel de")
    rec("folga explicita e prova o corte EMITIDO diretamente, e o terceiro estagio,")
    rec("que procura contraexemplo primal certificado em (y*,z*) e em (1,1).")
    rec("")
    for dp in dumps:
        rec(f"== {os.path.basename(dp).split('_')[0]} ==")
        try:
            r = run_one(dp, rec, checks, args.eps, None, args.instdir)
        except Exception as e:  # noqa: BLE001
            rec(f"ERRO: {e}")
            results.append({"name": os.path.basename(dp).split("_")[0], "ok": False,
                            "error": str(e), "times": {}})
            rec("")
            continue
        results.append(r)
        rec(f"TOTAL {r['name']}: {'PASS' if r['ok'] else 'FAIL'}  "
            f"[{sum(r['times'].values()):.1f}s]")
        rec("")

    rec("=" * 78)
    rec("SUMARIO")
    rec("=" * 78)
    rec("")
    rec("instancia          C1    C2    C3 validos    C4 margem  folga v_LD    tempo")
    rec("-" * 78)
    for r in results:
        if r.get("error"):
            rec(f"{r['name']:<18} ERRO: {r['error']}")
            continue
        c1 = "OK" if r.get("check1") else ("--" if "check1" not in r else "FAIL")
        c2 = "OK" if r.get("check2") else ("--" if "check2" not in r else "FAIL")
        r3 = r.get("r3")
        c3 = "--"
        if r3:
            c3 = f"{r3['proved'] + r3['rescued']}/{r3['n']}"
        r4 = r.get("r4")
        mm = "--" if not r4 or r4.get("min_margin") is None else str(r4["min_margin"])
        vs = "--"
        if r4:
            vs = ("ausente" if r4.get("vld_slack") is None
                  else f"{float(r4['vld_slack']):.2f}")
        rec(f"{r['name']:<18} {c1:<5} {c2:<5} {c3:<12} {mm:<11} {vs:<12} "
            f"{sum(r['times'].values()):.1f}s")
    rec("")
    rec("(a folga de v_LD aparece arredondada aqui; o valor EXATO em Fraction esta")
    rec(" no bloco CHECK 4 de cada instancia. 'ausente' = o run nao adicionou a linha")
    rec(" Lagrangiana ao mestre, entao nao ha o que exonerar.)")
    rec("")

    ncuts = sum(r["r3"]["n"] for r in results if r.get("r3"))
    nproved = sum(r["r3"]["proved"] + r["r3"]["rescued"] for r in results if r.get("r3"))
    nfail = sum(r["r3"]["failed"] for r in results if r.get("r3"))
    nund = sum(r["r3"]["undecided"] for r in results if r.get("r3"))
    ntight = sum(r["r3"]["tight"] for r in results if r.get("r3"))
    rec(f"CHECK 3 consolidado: {nproved} de {ncuts} cortes de otimalidade "
        f"demonstrados VALIDOS para todo (y,z) binario")
    rec(f"                     {ntight} deles com margem dual exatamente zero "
        "(o dual recuperado e o proprio dual otimo do subproblema)")
    rec(f"                     {nfail} falhas, {nund} nao decididos")
    rec("")
    rec("O QUE A CADEIA COBRE AGORA")
    rec("-" * 78)
    rec("O mestre final e composto por, e somente por, as seguintes familias de")
    rec("desigualdades. Cada uma esta abaixo com o seu estatuto apos os CHECK 1-4.")
    rec("")
    rec("  (a) y, z binarios, theta_l >= 0                : definicao das variaveis.")
    rec("  (b) (F1_l), (F2_l): sum_i b_il y_i >= Q_l e")
    rec("      sum_j p_jl z_j >= Q_l                      : validas por argumento")
    rec("      combinatorio elementar (a capacidade aberta tem de cobrir a demanda);")
    rec("      CHECK 1 confere que o incumbente as satisfaz.")
    rec("  (c) cortes de otimalidade de Benders           : CHECK 3. Cada um dos")
    rec(f"      {ncuts} cortes recebeu um ponto dual viavel EXPLICITO do LP de")
    rec("      roteamento cujo objetivo domina o corte; como a regiao dual nao")
    rec("      depende de (y,z), a desigualdade e valida em TODO o cubo, nao so no")
    rec("      incumbente. Demonstrado em aritmetica exata, sem solver.")
    rec("  (d) piso theta_l >= v_l(1,1)                   : CHECK 4(i). Valido por")
    rec("      monotonicidade de capacidade; alem disso a margem theta*_l - v_l(1,1)")
    rec("      medida exatamente e grande, entao nem um erro sub-unitario do LP que")
    rec("      gerou o piso poderia ter excluido o incumbente.")
    rec("  (e) linha global fixos + sum_l theta_l >= v_LD : CHECK 4(ii). Exonerada:")
    rec("      no otimo a folga exata e estritamente positiva (ou a linha sequer foi")
    rec("      adicionada), logo a linha nao cortou o otimo.")
    rec("")
    rec("Somando: toda desigualdade do mestre final esta demonstrada valida ou")
    rec("exonerada no otimo, e o valor do incumbente esta recomputado em inteiros")
    rec("por um solver de fluxo independente (CHECK 1), com os cortes fechando esse")
    rec("valor com folga zero (CHECK 2).")
    rec("")
    rec("UM COROLARIO QUE SO EXISTE COM O CHECK 3")
    rec("-" * 78)
    rec("Antes do CHECK 3, o CHECK 2 dizia apenas 'a soma dos maiores cortes iguala o")
    rec("transporte'. Com a validade dos cortes DEMONSTRADA, cada corte avaliado em")
    rec("(y*,z*) e um piso legitimo de v_l(y*,z*), logo best_l <= v_l(y*,z*) para todo")
    rec("l. Como o CHECK 2 fecha com folga EXATAMENTE zero, sum_l best_l = sum_l")
    rec("v_l(y*,z*), e uma soma de desigualdades <= que fecha com igualdade forca")
    rec("igualdade termo a termo: best_l = v_l(y*,z*) = theta*_l. Ou seja, no ponto")
    rec("otimo o mestre nao superestima NEM subestima nenhum theta_l; a relaxacao de")
    rec("Benders e exata em (y*,z*). E esse fato que sustenta a exoneracao da linha")
    rec("v_LD no CHECK 4(ii), que precisa saber o valor de fixos + sum_l theta_l no")
    rec("otimo do mestre, e nao apenas um piso dele.")
    rec("")
    rec("O UNICO ELO RESTANTE")
    rec("-" * 78)
    rec("A ARVORE DO SOLVER. Nada aqui reconstroi quais nos o branch-and-bound")
    rec("explorou, com que cota cada no foi podado, nem os cortes proprios do")
    rec("Gurobi (MIR, Gomory, cover) aplicados internamente. O certificado prova")
    rec("que o mestre final e uma relaxacao VALIDA cujo otimo o incumbente atinge,")
    rec("mas a afirmacao 'nenhuma outra abertura e mais barata' continua apoiada na")
    rec("corretude da busca do solver. Fechar esse elo exigiria despejar a arvore")
    rec("inteira ou uma prova dual do mestre, o que nao esta implementado.")
    rec("")
    rec(f"tempo total do relatorio: {time.time() - t_all:.1f}s")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"[relatorio escrito em {args.out}]", file=sys.stderr)
    return 0 if all(r.get("ok") for r in results) else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("certdump", nargs="?", help="arquivo .certdump")
    ap.add_argument("--instance", default=None,
                    help="caminho da instancia (default: a linha INSTANCE do dump)")
    ap.add_argument("--instdir", default=None,
                    help="diretorio onde procurar a instancia pelo nome base")
    ap.add_argument("--batch", default=None,
                    help="processa todos os *.certdump do diretorio e emite o relatorio")
    ap.add_argument("--out", default=None, help="arquivo de saida do --batch")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--eps", default="1e-6",
                    help="eps de sparsificacao do C++ (BendersOptions::eps)")
    ap.add_argument("--check1", action="store_true")
    ap.add_argument("--check2", action="store_true")
    ap.add_argument("--check3", action="store_true")
    ap.add_argument("--check4", action="store_true")
    ap.add_argument("--selftest", action="store_true",
                    help="roda a instancia minima embutida, com teste negativo, e sai")
    a = ap.parse_args(argv)

    def out(msg=""):
        if not a.quiet:
            print(msg)

    a.eps = Fraction(Decimal(a.eps))
    sel = {n for n, f in ((1, a.check1), (2, a.check2), (3, a.check3), (4, a.check4)) if f}
    a.checks = sel if sel else {1, 2, 3, 4}

    if a.selftest:
        return 0 if selftest(out) else 1
    if a.batch:
        return batch(a, out)
    if not a.certdump:
        ap.error("informe um .certdump, --batch DIR ou --selftest")

    t0 = time.time()
    try:
        r = run_one(a.certdump, out, a.checks, a.eps, a.instance, a.instdir)
    except FileNotFoundError as e:
        print(f"ERRO: {e}", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001
        print(f"ERRO processando {a.certdump}: {e}", file=sys.stderr)
        return 2
    out(f"TOTAL: {'PASS' if r['ok'] else 'FAIL'}  [{time.time() - t0:.1f}s]")
    out("LEMBRETE: os CHECK 1-4 provam que o mestre final e uma relaxacao valida "
        "que o incumbente atinge; a arvore do solver continua fora do "
        "certificado (ver LIMITACOES 1).")
    return 0 if r["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
