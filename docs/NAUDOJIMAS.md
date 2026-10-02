# Naudojimas

Šis dokumentas paaiškina, **ką rodo įrankis ir kaip tai skaityti**.
Diegimas aprašytas atskirai — žr. DIEGIMAS.md.

---

## Svarbiausia iš karto: ką šis įrankis daro ir ko nedaro

Įrankis **skaito jau paruoštą kandidatų rinkinį** (`.bcf` arba `.vcf` failą) ir
parodo, kokie duomenys BAM faile tuos kandidatus pagrindžia arba nepagrindžia.

**Jis pats kandidatų neieško.** Kandidatų rinkinį pagamina atskira programa —
`delly` — ir tai yra visiškai atskiras žingsnis, kurio čia nėra ir kuris
**nebuvo įdiegtas** kartu su šiuo įrankiu.

> Jeigu skaitote tik šį dokumentą, nesusidarykite įspūdžio, kad įdiegėte visą
> analizės grandinę. Įdiegėte **antrąją jos dalį** — skaitytuvą.

### Ką kainuoja pirmoji dalis (kandidatų gaminimas)

`delly` reikia dalykų, kurių šiam įrankiui nereikia: **referencinio genomo**
(apie 3 GB), pačios `delly` programos ir nemažai skaičiavimo laiko.

Pavyzdinė komanda (taip buvo pagamintas demonstracinis rinkinys):

```bash
delly sr \
  -g GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -x human.hg38.excl.tsv \
  -o meginys.bcf \
  --threads 1 \
  meginys.bam
```

Kiek tai trunka — išmatuota:

| BAM dydis | Trukmė | Didžiausias atminties poreikis |
|---|---|---|
| ~1,6 GB (dvi chromosomos) | **39 s** | mažas |
| ~40 GB (visas genomas) | **46 min** ir **2 val. 4 min** | apie **1,2 GB** |

> *Papildyta 2026-09-26.* 46 min ir 2 val. 4 min — pirmojo paleidimo laikai. Tuos pačius
> mėginius 2026-09-25 apdorojus iš naujo ta pačia `delly` versija ir viena gija, tai užtruko
> 1 val. 2 min ir 1 val. 45 min, o atminties prireikė 1 158 ir 1 187 MiB
> (`stage1_igv_assistant/results/patient_rerun_2026-09.json`): trukmė priklauso nuo
> kompiuterio apkrovos, atminties poreikis beveik nesikeičia.

Tai daroma **vieną kartą** kiekvienam mėginiui. Gautą `.bcf` failą (paprastai
apie 130 KB) po to galima skaityti šiuo įrankiu kiek nori kartų.

---

## Darbo eiga

Naršyklėje atverkite **http://127.0.0.1:8765**. Puslapyje trys žingsniai, matomi
viršuje kaip skirtukai: **1 Candidates**, **2 Evidence** ir **3 Ask the assistant**.
Viršutinėje juostoje — mėginio pasirinkimas **Sample** ir du mygtukai: **Limits**
(ko įrankis pasakyti negali) ir **Call log** (visi šios sesijos įrankių kvietimai).
Šalia skirtukų — būsena: kiek įrankių patikrinta paleidžiant, kiek prieinama
vietinių ir debesijos modelių, ar rastas IGV.

*Pastaba, 2026-10-02:* puslapis perdarytas (Phase 24). Ankstesnis puslapis su tais
pačiais įrankiais tebėra adresu **http://127.0.0.1:8765/classic**.

### Testiniai ir privatūs duomenys

Prie mėginio pavadinimo matomas ženkliukas **Test data** arba **Private data**.

| Ženkliukas | Kas taip žymima |
|---|---|
| **Test data** | failai, kuriuos įrankis rado pats viešų duomenų kataloge (`data_dir`), ir etiketės, kurias konfigūracijos failas skiltyje `[test_data]` paskelbė testiniais duomenimis |
| **Private data** | visa, kas užregistruota aiškiai: komandinėje eilutėje (`--dataset`, `--candidates`) arba konfigūracijos failo skiltyse `[datasets]` ir `[candidates]` |

Puslapis visada atsidaro su testiniu mėginiu. Privatus mėginys savaime
neįkeliamas — tik jį pasirinkus. Ir tada rodomi **tik suvestiniai skaičiai**:
kiek jungčių lieka po kiekvieno filtro. Kandidatų sąrašas su pozicijomis
atsiranda tik paspaudus **Show the … candidates** (iki tol puslapis iš serverio
pozicijų net negauna). Taip privataus mėginio pozicijos nepamatomos atsitiktinai,
pvz., sąraše **Sample** paspaudus rodyklės klavišą. Kada privatūs duomenys gali
palikti kompiuterį, aprašyta skyriuje „Asistentas“.

Viena išimtis — palyginimas su kitu mėginiu (**Mark junctions also found in**).
Jam serveris sudaro abiejų mėginių filtruotus sąrašus; puslapis ir tada gauna tik
skaičius, bet tie du kvietimai su pozicijomis lieka žurnale **Call log**. Be to,
jei palyginimui pasirinktas privatus mėginys, kito mėginio eilutės su žyme
„also in …“ parodo, kurias jungtis tas privatus mėginys turi bendras — taigi ir
jo pozicijas (laukelyje **within … bp** nurodyto atstumo tikslumu), nors jo
paties sąrašas neatvertas.

*Pastaba, 2026-10-02 (Phase 25):* paslėptas privataus mėginio sąrašas, leidimo
apimtis ir platesnis klausimo patikrinimas (žr. „Asistentas“) pridėti po pirmojo
bandymo su tikrais duomenimis.

### 1. Candidates — kandidatų sąrašas

Sąraše **Sample** pasirinkite mėginį. Įrankis įkelia jo kandidatų rinkinį (sąraše
matomos **etiketės**, ne failų keliai) ir iš karto pritaiko filtrus. Pagal
nutylėjimą rodomos jungtys **tarp dviejų chromosomų** (BND). Pakeitus bet kurį
filtrą, sąrašas perskaičiuojamas iškart — mygtuko spausti nereikia.

| Filtras | Ką daro |
|---|---|
| **Type** | palieka vieno tipo jungtis: tarp dviejų chromosomų (BND), delecijas (DEL), duplikacijas (DUP), inversijas (INV) — arba visas |
| **Caller marked PASS** | tik tas, kurias `delly` pažymėjo PASS |
| **Read pairs at least** | mažiausias porinių skaitinių skaičius (PE), kurį pranešė `delly` |
| **Split reads at least** | mažiausias perskeltų skaitinių skaičius (SR), kurį pranešė `delly` |
| **Main chromosomes only** | abu jungties galai pagrindinėse chromosomose |
| **Skip known problem regions** | atmeta jungtis, kurių galas patenka į delly neįtraukiamų sričių šabloną |
| **Mark junctions also found in … within … bp** | pažymi jungtis, kurios yra ir kitame mėginyje (abu galai ne toliau nei nurodyta) |
| **Only those that touch a gene** | palieka tik persitvarkymus, kurių lūžis yra gene arba perkeltame segmente yra genas (reikia genų lentelės) |

**Visi kandidatai įvertinami iš karto.** Įrankis sujungia jungtis, kurios priklauso
tam pačiam persitvarkymui (pvz., abi subalansuotos translokacijos jungtys), kiekvienam
suskaičiuoja palaikančius skaitinius tiesiai iš BAM failo, įvardija genus lūžių
vietose ir surikiuoja sąrašą **pagal palaikančių skaitinių skaičių**. Balų nėra.

Lentelės stulpeliai:

| Stulpelis | Kas jame |
|---|---|
| **Rearrangement** | pavadinimas ISCN principu, pvz., `t(19;22)(q13.33;q12.2)` (juostos apskaičiuotos iš koordinačių pagal GRCh38), ir pobūdis: **Both junctions found** (rastos abi subalansuotos translokacijos jungtys), **One junction**, **Segment moved** (perkeltas segmentas), **Complex** |
| **Breakpoints** | lūžių koordinatės abiejose chromosomose |
| **Genes at the breakpoints** | genai lūžių vietose (kuris intronas ar egzonas), **OMIM** žymė, jei genas susijęs su liga; jei geno nėra — artimiausias |
| **Supporting reads** | kiek skirtingų skaitinių palaiko persitvarkymą: porinių skaitinių, kurių vienas galas prie vieno lūžio, kitas — prie kito, ir perskeltų skaitinių, kurių dalys prie abiejų lūžių; abu galai patikimai nusėdę (MAPQ ≥ 20) |
| **Read with care** | perspėjimai: lūžis centromeroje ar heterochromatine, daug skaitinių rodo į kitas vietas, daug dviprasmiškai nusėdusių skaitinių, jungtis yra ir kitame mėginyje, abipusė jungtis faile yra, bet atmesta filtrų |

Eilutę paspaudus (arba **Review**) atidaromas 2 žingsnis. Po lentele —
**The caller's … junctions**: tos pačios jungtys taip, kaip jas pranešė `delly`.

**How the filters narrow the list** (po sąrašu) rodo, kiek jungčių lieka po
kiekvieno žingsnio, iš kur kiekvieno žingsnio riba (**set by the tool**,
**author's choice**, **reference file**) ir kiek jis pašalino grandinėje ir vienas
(**(−23 alone)**). Kai **Main chromosomes only** ar **Skip known problem regions**
nieko nepašalina, puslapis žodžiais paaiškina kodėl: faile tokių jungčių nėra arba
jas jau pašalino ankstesni žingsniai. `delly`, paleistas su tuo pačiu neįtraukiamų
sričių šablonu, tokias jungtis praleidžia pats.

### 2. Evidence — vienas persitvarkymas

Svarbiausia dalis — **schema** (**What it does to the chromosomes**):

- viršuje — abi chromosomos su juostomis ir centromera, lūžio vieta (koordinatė,
  juosta, genas);
- spalvotu kontūru pažymėtas **gabalas, kuris išvyksta** iš kiekvienos chromosomos,
  o rodyklė rodo, **prie kurios chromosomos ir kurioje vietoje jis prisijungia**;
- apačioje (**After the rearrangement**) — chromosomos, kurias sukuria kiekviena
  jungtis (pvz., `der(19)`), su kiekvieno gabalo pradžios ir pabaigos
  koordinatėmis.

Toliau:

- **Genes** — genai kiekviename lūžyje (egzonas ar intronas pagal kanoninį
  transkriptą, grandinė, OMIM numeriai su nuorodomis), genai tarp dviejų gretimų
  jungčių ir kiek genų (ir kiek su liga susijusių) yra išvykstančiame gabale;
  **List them** juos išvardija. Genai skaitomi iš lentelės šiame kompiuteryje —
  pozicijos niekur nesiunčiamos.
- **The reads at the two ends** — dvi dėžutės (kaip du IGV langai): kairėje
  viena chromosoma, dešinėje kita. Rodomi **tik nenormalūs skaitiniai**; normalūs
  paslėpti (jų skaičius parašytas). Spalva rodo, ką skaitinys sako apie jungtį:
  viena spalva — skaitiniai, jungiantys abu galus taip, kaip viena jungtis, kita —
  kaip kita jungtis; pilki — mate kitoje chromosomoje; geltoni — mate nenusėdęs;
  žali — neteisinga porų orientacija; raudoni — per didelis atstumas tarp porų;
  tuščiaviduriai (balti) — dviprasmiškai nusėdę (MAPQ < 20), jie nupiešti, bet
  neskaičiuojami. Punktyrinė linija — lūžis; viršuje genai. Užvedus pelę ant
  skaitinio matyti jo vardas, vieta, CIGAR ir kur jo kita dalis.
- **Reads for …** — patys skaitiniai lentelėmis (vardas, pozicija, grandinė,
  CIGAR, kokybė abiejuose galuose).
- **What the caller reported** — `delly` įrašai ir jų PE/SR palyginti su čia
  suskaičiuotais skaitiniais.
- **IGV images at …** — IGV paveikslėliai (žr. žemiau); nebūtini.

**Bet kurią kitą vietą** galima patikrinti dešinėje viršuje: pasirinkite, kurio
mėginio skaitinius naudoti (**Reads from**), įrašykite vietą (pvz.,
`chr20:33,700,000`) ir paspauskite **Check position**. Puslapis parodo, į kurias
chromosomas rodo tos vietos skaitiniai, kur jie nukirpti ir kokie genai ten yra.
Perskelti skaitiniai duoda tikslias partnerio pozicijas: mygtukas **Junction with …**
atidaro tą vietą kaip jungtį su abiem galais. Puslapis aiškiai parašo **Typed in by
hand**: jos nepatvirtino joks kandidatų rinkinys.

### Genų lentelė (vieną kartą)

Genai rodomi, kai sukurta vietinė genų lentelė:

```bash
mkdir -p ~/public_data/annotation && cd ~/public_data/annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz
python3 ~/rare-disease-diagnosis-assistant/scripts/make_gene_table.py gencode.v46.basic.annotation.gtf.gz genes_grch38.tsv.gz
wget https://omim.org/static/omim/data/mim2gene.txt                         # nebūtina: OMIM genų numeriai
wget -O genes_to_disease.txt https://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt   # nebūtina: genas → liga
```

Numatytosios vietos: `~/public_data/annotation/genes_grch38.tsv.gz`,
`mim2gene.txt`, `genes_to_disease.txt` (arba konfigūracijos failo `[paths]`
raktai `gene_table`, `omim_mim2gene`, `gene_disorders`). Paleidimo eilutėse
matyti, ar lentelė rasta. Be lentelės įrankis veikia, tik genų neįvardija.

### IGV paveikslėliai

**IGV images at …** paleidžia IGV, kuris nupiešia po paveikslėlį kiekvienam
matavimui. Tai trunka **kelias minutes**. IGV rodo **visus** skaitinius, ir
normalius; aukščiau esančios dvi dėžutės rodo tik nenormalius. Paveikslėlius mato
tik žmogus — asistentas jų nemato. Privatiems duomenims IGV prašo `igv.org` tos
srities referencinės sekos: srities koordinatės palieka kompiuterį, skaitiniai —
ne.

### 3. Ask the assistant

Žr. skyrių „Asistentas“. Į klausimą „Which rearrangements in SAMPLE have the most
read support, and which genes do they break?“ asistentui nurodyta atsakyti vienu
įrankio kvietimu (`review_candidates`), kuris įvertina visus kandidatus.

*Pastaba, 2026-10-03 (Phase 26):* po susitikimo su vadove balų sistema iš puslapio
pašalinta; skyriai apie balus žemiau aprašo ankstesnį puslapį
(**http://127.0.0.1:8765/classic**).

---

## Keturi įrodymų sluoksniai

*Šis ir du kiti skyriai su balais aprašo ankstesnį puslapį (`/classic`). Naujame
puslapyje balų nėra.*

| Sluoksnis | Ką skaičiuoja | Ką reiškia |
|---|---|---|
| **Poriniai skaitiniai** (discordant pairs) | poras, kurių antras galas nusėdo kitoje chromosomoje | rodo, su kuo sujungta |
| **Nukirsti skaitiniai** (soft-clipped) | skaitinius, kurių galas „nukerpamas“ toje pačioje vietoje | rodo **tikslią** lūžio vietą |
| **Perskelti skaitiniai** (split reads) | skaitinius, kurių dalis nusėdo kitur | pats tiesiausias įrodymas |
| **Skaitymo gylis** (read depth) | ar toje vietoje sumažėjo skaitymo gylis (taškai skiriami tik už kritimą) | rodo iškritas; padaugėjimas (dublikacija) taškų negauna |

Kiekvienas vertinamas atskirai nuo 0 iki 25, iš viso 0–100.
Vertinimas: **70 ir daugiau — „strong“**, **40 ir daugiau — „moderate“**,
daugiau nei 0 — „weak“.

Prie kiekvieno sluoksnio skaičiaus ir prie bendro įverčio yra nuoroda **source**,
atidaranti **tikslų įrankio atsakymą**, kuris tą skaičių grąžino. Pasiekiama riba
nėra įrankio atsakymas: sąsaja ją išveda iš vertinimo pakopų ir stebėtų reikšmių
(dalių ir nukirptų skaitinių skaičiaus), skaičiuodama tuos pačius sluoksnius kaip
įvertis. Visi sesijos kvietimai — mygtuku **Call log**.

---

## „QUALITY-LIMITED“ — tai ne blogas įvertinimas

Jei per daug skaitinių toje vietoje nusėdo nepatikimai (daugiau nei 40 %
žemiau MAPQ 20), įrankis **neskaičiuoja** bendro įvertinimo ir parašo
`QUALITY-LIMITED`.

> **Tai reiškia „nepamatuota“, o ne „pamatuota ir mažai“.**

Skirtumas esminis. „weak“ reiškia: pažiūrėjome ir beveik nieko neradome.
`QUALITY-LIMITED` reiškia: **negalėjome pažiūrėti** — ši genomo vieta tokia
pasikartojanti, kad nežinia, ar skaitiniai apskritai iš čia.

Puslapyje vietoje skaičiaus rodoma **Withheld**. Keturi atskiri sluoksniai vis
tiek rodomi. Skaitykite juos.

Jei nė vieno matavimo toje vietoje padaryti nepavyko (pvz., ten nėra skaitinių),
rodoma **No score** ir „not assessable“, o lentelė kiekvienam matavimui
paaiškina kodėl.

---

## Kodėl subalansuota translokacija čia negauna „strong“

Tai aritmetika, ne duomenų trūkumas.

Subalansuotoje translokacijoje **medžiagos nei pridedama, nei atimama**.
Todėl skaitymo gylio sluoksnis teisingai duoda **0 taškų iš 25** — ir tai
teisingas atsakymas, ne klaida.

Antra: jei pertvarkyta **tik viena iš dviejų** chromosomos kopijų
(heterozigotinis atvejis), maždaug **pusė** skaitinių toje vietoje ateina iš
**sveikosios kopijos**. Todėl porinių skaitinių dalis niekada nepasiekia 0,5,
kurios reikalauja aukščiausia to sluoksnio pakopa.

Sudėjus: du sluoksniai iš keturių yra apriboti iš anksto, ir aukščiausia
juosta pasidaro **nepasiekiama** — nebent gylio matavimas klaidingai duotų taškų.

Įrankis tai apskaičiuoja **kiekvienai vietai atskirai** ir parodo laukuose
`attainable_here` (didžiausias čia pasiekiamas įvertinimas) ir `strong_band`
(nuo kiek prasideda „strong“). Puslapyje tai brūkšninis ženklas **reachable here**
ant įverčio juostos. Jei pirmasis mažesnis už antrąjį — „strong“ čia
nepasiekiamas, kad ir kokie geri būtų skaitinių duomenys.

Delecijoms ir duplikacijoms šis ženklas nerodomas: jis remiasi prielaida, kad
medžiagos nepridėta ir neatimta, o delecija ir duplikacija skaitymo gylį keičia.
Ranka įvestai vietai jis rodomas su sąlyga „jei tai subalansuotas pertvarkymas“.

*Pastaba, 2026-09-30:* nuo pataisymo 25d12bb riba skaičiuojama tik pagal tuos sluoksnius, kuriuos skaičiuoja pats įvertis (tinkami sluoksniai, atėmus neįvertinamus). Jei nurodoma mažiau nei keturi sluoksniai, atsakyme atsiranda laukas `ceiling_counted_layers`, o paaiškinimas juos įvardija. Kai skaičiuojami visi keturi, laukai tokie patys kaip anksčiau.

> Subalansuotą translokaciją vertinkite pagal **keturis atskirus matavimus**,
> o ne pagal juostą, į kurią pateko bendras skaičius.

---

## Asistentas

Skirtukas **Ask the assistant**. Asistentas — kalbos modelis, kuris duomenis
pasiekia **tik per tuos pačius įrankius**, kuriais naudojasi puslapis. Kiekvienas
skaičius jo galutiniame atsakyme patikrinamas su tuo, ką tie įrankiai grąžino:
patvirtinti skaičiai pabraukti žaliai, o skaičius, kurio negrąžino joks įrankis,
pažymimas raudonai. Tikrinami tik skaičiai, ne žodžiai: tokius žodžius kaip
„stipriausias“ modelis renkasi pats. Todėl skaitykite ir įrankių rezultatus.

Asistentui siūlomi tie patys įrankiai kaip puslapiui, **išskyrus bendrą balą**
(`breakpoint_evidence_summary`): įrodymus jis aprašo skaitinių skaičiais. Į klausimą
„kurie kandidatai stipriausi“ jam nurodyta atsakyti vienu kvietimu `review_candidates`, kuris
įvertina visus kandidatus ir grąžina juos surikiuotus. Genus jis įvardija įrankiu
`genes_near`, kuris skaito vietinę genų lentelę; kai lentelė sukurta, internetinė
genų paieška (Ensembl) jam nebesiūloma.

1. **Model** — pasirinkite modelį:
   - **On this computer (Ollama)** — vietiniai modeliai (pvz., `qwen2.5:7b`).
     Klausimas ir duomenys kompiuterio nepalieka (išskyrus IGV ir, jei genų lentelė
     nesukurta, internetinę genų paiešką).
   - **Cloud (Anthropic)** — **Claude Sonnet 5** ir **Claude Opus 5**. Rodomi tik
     tada, kai yra API raktas. Klausimas ir kiekvienas įrankio atsakymas, kurio
     modelis paprašo, siunčiami į Anthropic serverius; atsakymo apačioje nurodyta
     jo kaina.
2. **Your question** — įrašykite klausimą arba paspauskite vieną iš
   **Suggested questions** (jie sudaromi pagal paskutinę atvertą vietą).
3. **Ask** (arba Ctrl+Enter). Vienu metu — vienas klausimas.

Atsakymo kortelėje: atsakymas; patikrinimo eilutė (kiek skaičių patikrinta arba
kiek nepatvirtinta); kiekvienas modelio padarytas įrankio kvietimas — viena eilutė
su santrauka, **Raw result** (visas įrankio atsakymas) ir **source** (įrašas
kvietimų žurnale); apačioje — trukmė, modelio žingsnių skaičius ir debesijos
modeliui kaina.

### Privatūs duomenys ir debesijos modeliai

Debesijos modeliui privatūs duomenys siunčiami **tik jūsų patvirtinimu ir tik
tam vienam klausimui**:

- Jei klausimas **mini privatų mėginį** arba **turi vietą, perskaitytą iš privačių
  duomenų, ar šalia jos** (iki 1 000 bazių; atpažįstami užrašai `33700000`,
  `33,700,000`, `33 700 000`, `33.700.000`, `33.7 Mb`, `33,7 Mb`, `33700 kb`), po
  klausimu atsiranda oranžinis langelis, o **Ask** neaktyvus, kol nepažymite
  **I have permission to send this private data to Anthropic**. Išsiuntus klausimą
  žymė nuimama: kitam klausimui leidimą reikia pažymėti iš naujo.
  Tikrinamas tik skaičius, ne chromosoma, todėl, įkėlus didelį privatų kandidatų
  rinkinį, langelis gali atsirasti ir dėl su juo nesusijusio skaičiaus.
  Neatpažįstama: skaičius, prie kurio be tarpo prirašytos raidės (`33700000bp`),
  `33.7M`, grupuota kb reikšmė (`33,700 kb`), mokslinis užrašas (`3.37e7`).
- Nepažymėjus modelis privačių etikečių **net nemato**: jam pasiūlomos tik
  testinių duomenų etiketės, o įrankio kvietimas su privačia etikete arba su iš
  privataus failo įkelto rinkinio identifikatoriumi atmetamas dar prieš jį
  vykdant (atsakyme — **Blocked to protect private data**). Serveris klausimo
  tekstą tikrina dar kartą, todėl, net ir apėjus puslapį, toks klausimas
  neišsiunčiamas.
- Failo kelias vietoje etiketės (taip pat ir neįtraukiamų sričių failo kelias)
  atmetamas **bet kuriam** modeliui, ir vietiniam: duomenys pasiekiami tik per
  etiketes.
- Leidimas apima **tik tuos privačius mėginius, kuriuos klausimas mini arba
  kurių pozicijas jis turi** — tuos pačius, kuriuos įvardija oranžinis langelis.
  Kitų privačių mėginių modelis nemato, o jų kvietimai atmetami.
- Pažymėjus atsakyme parašoma, kas išsiųsta ir ką leidimas apėmė: **Sent to
  Anthropic with your permission: tool results for …**.

Siųskite privačius duomenis tik turėdami duomenų savininko leidimą. Vietiniams
modeliams šie apribojimai, išskyrus failo kelio atmetimą, netaikomi — duomenys
lieka kompiuteryje.

**Genų paieška** (`gene_at_locus`) siunčia chromosomą ir poziciją į Ensembl
(`rest.ensembl.org`) — ir naudojant vietinį modelį. Tokie kvietimai pažymėti
**uses internet**. Privačioms vietoms ją naudokite tik tada, kai tai leidžiama.

---

## Demonstracinio rinkinio ribos

Demonstracinis rinkinys yra **1,6 MB**. Jame yra du mėginiai, o jų BAM failai
**iškarpyti**: palikti tik reikalingi ruožai, o ne visos chromosomos. Taip
daroma sąmoningai — pilni failai svertų po 1,6 GB kiekvienas.

| Mėginys | Kas jame palikta | Ką parodo |
|---|---|---|
| `DEMO_CLEAN` | `chr20:200 000` ir `chr21:14 100 000`, po ±10 000 bazių | trys iš keturių sluoksnių; įvertis 40,0/100 „moderate“; aukščiausia juosta čia nepasiekiama (daugiausia 57,5, reikia 70) |
| `DEMO_REPEAT` | `chr20:25 800 000` ir `chr21:7 600 000`, po ±10 000 bazių | `QUALITY-LIMITED` — įvertis sulaikytas |

> *Pataisyta 2026-09-25.* Čia buvo parašyta „visi keturi sluoksniai; įvertis
> 47,5/100“ — tai pirmojo implantų rinkinio reikšmė; to rinkinio įrašai prarasti.
> Rinkinys perkurtas, ir dabar `make_demo_bundle.py` šias reikšmes kiekvieną kartą
> išmatuoja pačiame iškarpytame faile ir įrašo į `DEMO.md`.

| Veikia visiškai | Veikia tik ruožuose |
|---|---|
| filtrų grandinė (896 sujungimai po sulydymo) | keturi įrodymų sluoksniai |
| dviejų rinkinių palyginimas | ranka įvesta koordinatė |

Filtrų grandinė ir palyginimas veikia **pilnai**, nes jie skaito tik `.bcf`
failus, o tie nesukarpyti.

Įvedus koordinatę už ruožų ribų, sluoksniai parodys, kad **duomenų ten nėra**
(`assessable: false`) — tai ne klaida, o sąžiningas atsakymas. Tikruose
duomenyse tokio apribojimo nėra.

> **Patikrinta:** tuose ruožuose iškarpytas failas duoda **tą patį** įvertį kaip
> ir pilnas — 40,0/100 (pirmajame implantų rinkinyje buvo 47,5/100). Tai nėra savaime suprantama: pirmasis bandymas davė
> 43,3, nes sluoksnių tinkamumas nustatomas pagal failo pradžios skaitinius, o
> iškarpytame faile jie yra kitokie. Tai ištaisyta.

---

## Ko šis įrankis pasakyti negali

**Jis tikrina pozicijas, bet jų neieško.** Kiekviena čia tikrinama vieta arba
atėjo iš kandidatų rinkinio, arba buvo įvesta ranka. Tikras lūžis, kurio
neaptiko `delly` ir kurio niekas neįvedė, čia **niekada nepasirodys**.

**Kontroliuotame teste rasti 16 iš 24 žinomų lūžio taškų.** Švarioje,
vienareikšmiškai skaitomoje sekoje — 8 iš 8. Šalia pasikartojančių sričių —
8 iš 8. Ten, kur seka skaitoma dviprasmiškai — **0 iš 8**. Visi praleisti buvo
prarasti kandidatų paieškos žingsnyje, ne filtruose.

> *Pataisyta 2026-09-25.* Čia buvo parašyta „14 iš 24“ ir „šalia pasikartojančių
> sričių 6 iš 8“ — tai pirmojo bandymo rezultatas; jo įrašai prarasti. Testas
> perkurtas kaip naujas eksperimentas. Skirtumą lemia IMP06: tada praleistas,
> dabar aptiktas; perlyginus jo skaitinius be `.alt` failo jis vėl prarandamas,
> taigi skirtumą paaiškina ALT-aware lyginimas, o ne atsitiktinė imtis.

**14 iš 16 ribų yra autoriaus sprendimas, o ne kalibruotos vertės.** Jos
nepatikrintos su patvirtintų teigiamų ir neigiamų atvejų rinkiniu. Prie
kiekvienos ribos parašyta, iš kur ji. „Autoriaus sprendimas“ skaitykite taip:
protingas žmogus būtų pasirinkęs kitaip, ir rezultatas būtų kitoks.

**Skaitymo gylio matavimas prie tokio padengimo nepatikimas.** Jo riba nustatyta
su maždaug dešimt kartų didesnio padengimo duomenimis. Vertinkite jį kaip
silpną papildomą požymį, ne daugiau.

> *Pataisyta 2026-09-27.* Perskeltų skaitinių sluoksnis anksčiau skaičiavo papildomą
> sulygiavimą ant „decoy“ kontigo (`chrUn_..._decoy`) kaip jungties partnerį. Viešuose
> NA12878 duomenyse vienoje vietoje (chr21:10 770 078) dėl to įvertis buvo 72,5 „strong“
> vietoje 55,0 „moderate“. Dabar tokie įrašai partneriais nelaikomi, bet rodomi atskirai
> (`decoy_partners`, `decoy_only_reads`), kad kartografavimo dviprasmiškumas liktų matomas.
> Jokia riba nepakeista.

**Įvertinimas nėra tikimybė.** Tai suskaidytas aprašas, kas suveikė, o ne
apskaičiuota tikimybė, kad variantas tikras.

**Tai magistro darbo prototipas, ne klinikinis įrankis.** Nė vienas čia matomas
skaičius netinka klinikiniam sprendimui be nepriklausomo patvirtinimo.
