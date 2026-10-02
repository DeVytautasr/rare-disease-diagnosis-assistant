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
neįkeliamas — tik jį pasirinkus. Kada privatūs duomenys gali palikti kompiuterį,
aprašyta skyriuje „Asistentas“.

### 1. Candidates — kandidatų sąrašas

Sąraše **Sample** pasirinkite mėginį. Įrankis įkelia jo kandidatų rinkinį (sąraše
matomos **etiketės**, ne failų keliai) ir iš karto pritaiko filtrus. Pakeitus bet
kurį filtrą, sąrašas perskaičiuojamas iškart — mygtuko spausti nereikia.

| Filtras | Ką daro |
|---|---|
| **Type** | palieka vieno tipo jungtis: translokacijas (BND), delecijas (DEL), duplikacijas (DUP), inversijas (INV) — arba visas |
| **Caller marked PASS** | tik tas, kurias `delly` pažymėjo PASS |
| **Read pairs at least** | mažiausias porinių skaitinių skaičius (PE) |
| **Split reads at least** | mažiausias perskeltų skaitinių skaičius (SR) |
| **Main chromosomes only** | abu jungties galai pagrindinėse chromosomose |
| **Skip known problem regions** | atmeta jungtis, kurių galas patenka į delly neįtraukiamų sričių šabloną |
| **Drop junctions also found in … within … bp** | atmeta jungtis, kurios yra ir kitame pasirinktame mėginyje (abu galai ne toliau nei nurodytas bazių skaičius) |

**How the filters narrow the list** rodo, kiek jungčių lieka po kiekvieno žingsnio.
Prie žingsnio pažymėta, iš kur jo riba: **set by the tool** (nustatyta įrankio),
**author's choice** (autoriaus sprendimas) arba **reference file** (iš
referencinio failo). Dešinėje — kiek liko; po skaičiumi — kiek šis žingsnis
pašalino grandinėje (pvz., **−11**) ir, jei skiriasi, kiek būtų pašalinęs vienas
(pvz., **(−23 alone)**). **no change** reiškia, kad grandinėje jis nieko nepašalino.

Šie du skaičiai dažnai skiriasi, ir tai svarbu. Jei filtras grandinėje nieko
nepašalino, bet vienas būtų pašalinęs 300, jis **nėra nenaudingas** — tiesiog
ankstesni filtrai tuos įrašus jau buvo pašalinę. Eilė turi reikšmės: tie patys
filtrai kita tvarka duoda tuos pačius galutinius kandidatus, bet kitokius
tarpinius skaičius. Juostų ilgiai — logaritminiu masteliu.

**Jei žingsnis neatliktas**, virš juostų atsiras raudonas užrašas su paaiškinimu.
Dažniausia priežastis — nerastas delly neįtraukiamų sričių šablonas. Tada
grandinė veikia, tik be to vieno žingsnio, ir apie tai pasako aiškiai.

Žemiau — **… candidates to review**, likusios jungtys. Abu jungties galai rodomi
genomo tvarka (**One end**, **Other end**), toliau tipas, porinių ir perskeltų
skaitinių skaičius ir `delly` žymė. Kai pasirinktas kitas mėginys, jame rastos
jungtys paslepiamos; pažymėjus **Also show the … found in …**, jos rodomos
pilkai su žyme „also in …“. **Comparison details** — abiejų rinkinių palyginimo
skaičiai: kiek jungčių yra abiejuose ir kiek tik viename. Naudinga ir tada, kai
tas pats mėginys apdorotas skirtingais nustatymais — matyti, ką pakeitimas
realiai pridėjo arba atėmė.

### 2. Evidence — įrodymai vienoje jungtyje

Kandidato eilutėje paspauskite **Review**. Įrankis paleidžia keturis įrodymų
sluoksnius abiejuose jungties galuose ir atidaro skirtuką **Evidence**:

- viršuje — jungtis (pvz., `chr20:200,000 ↔ chr21:14,100,001`), ką apie ją
  pranešė `delly`, ir iš kurio mėginio skaitiniai;
- schema — kur abu galai yra chromosomose;
- po skydelį kiekvienam galui: **bendras įvertis** iš 100 su juosta, ženklas
  **reachable here** (žr. „Kodėl subalansuota translokacija čia negauna
  „strong““) ir lentelė — kiekvienas matavimas, jo reikšmė ir taškai;
- **Show IGV images** — IGV paveikslėliai (žr. žemiau);
- **Ask the assistant about this position** — pereina į 3 žingsnį su paruoštu
  klausimu apie šią vietą.

Prie bendro įverčio, prie kiekvieno matavimo ir prie kokybės eilutės yra nuoroda
**source**. Ji atidaro tikslų įrankio kvietimą: ką įrankis gavo ir ką grąžino.
Ženklas **reachable here** tokios nuorodos neturi: tai ne įrankio atsakymas, o
puslapio išvestas dydis (žr. „Keturi įrodymų sluoksniai“).

**Bet kurią kitą vietą** galima patikrinti dešinėje viršuje: pasirinkite, kurio
mėginio skaitinius naudoti (**Reads from**), įrašykite vietą (pvz.,
`chr20:33,700,000`) ir paspauskite **Check position**. Vieta **nebūtinai turi būti
iš kandidatų rinkinio**. Tada puslapis aiškiai parašo **Typed in by hand**: jos
nepatvirtino joks kandidatų rinkinys. Tai apsauga nuo savęs apgaudinėjimo —
radus „įrodymų“ ranka įvestoje vietoje, tai dar nereiškia, kad ten yra tikras
lūžis.

### IGV paveikslėliai

**Show IGV images** paleidžia IGV, kuris nupiešia po paveikslėlį kiekvienam
matavimui: porinių skaitinių, nukirptų skaitinių, perskeltų skaitinių ir
skaitymo gylio vaizdą. Paveikslėliai atsiranda per visą plotį po skydeliais;
paspaudus paveikslėlis padidinamas. Tai trunka **kelias minutes**, nes IGV
kiekvienam paveikslėliui iš interneto įkelia genomą ir genų takelį. Vienu metu
piešiamas vienas rinkinys: kol jis piešiamas, kiti IGV mygtukai neaktyvūs.
Paveikslėlius mato tik žmogus — asistentas jų nemato.

### 3. Ask the assistant

Žr. skyrių „Asistentas“.

---

## Keturi įrodymų sluoksniai

| Sluoksnis | Ką skaičiuoja | Ką reiškia |
|---|---|---|
| **Poriniai skaitiniai** (discordant pairs) | poras, kurių du galai nusėdo ne ten, kur turėtų | rodo, su kuo sujungta |
| **Nukirsti skaitiniai** (soft-clipped) | skaitinius, kurių galas „nukerpamas“ toje pačioje vietoje | rodo **tikslią** lūžio vietą |
| **Perskelti skaitiniai** (split reads) | skaitinius, kurių dalis nusėdo kitur | pats tiesiausias įrodymas |
| **Skaitymo gylis** (read depth) | ar padaugėjo/sumažėjo medžiagos | rodo iškritas ir dublikacijas |

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
pažymimas raudonai. Tikrinami tik skaičiai, ne žodžiai: modelis gali pavadinti
„strong“ tai, ką įrankis įvertino „weak“. Todėl skaitykite ir įrankių rezultatus.

1. **Model** — pasirinkite modelį:
   - **On this computer (Ollama)** — vietiniai modeliai (pvz., `qwen2.5:7b`).
     Klausimas ir duomenys kompiuterio nepalieka (išskyrus genų paiešką, žr. žemiau).
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
  duomenų** (pvz., ką tik patikrintą koordinatę), po klausimu atsiranda oranžinis
  langelis, o **Ask** neaktyvus, kol nepažymite **I have permission to send this
  private data to Anthropic**. Išsiuntus klausimą žymė nuimama: kitam klausimui
  leidimą reikia pažymėti iš naujo.
- Nepažymėjus modelis privačių etikečių **net nemato**: jam pasiūlomos tik
  testinių duomenų etiketės, o įrankio kvietimas su privačia etikete, su iš
  privataus failo įkelto rinkinio identifikatoriumi arba su tiesiogiai nurodytu
  failo keliu atmetamas dar prieš jį vykdant (atsakyme — **Blocked to protect
  private data**). Serveris klausimo
  tekstą tikrina dar kartą, todėl, net ir apėjus puslapį, toks klausimas
  neišsiunčiamas.
- Pažymėjus atsakyme parašoma, kas išsiųsta: **Sent to Anthropic with your
  permission: tool results for …**. Leidimas galioja visam klausimui, ne vienam
  mėginiui: tam klausimui modeliui pasiūlomos **visos** privačios etiketės, ne tik
  paminėtos klausime.

Siųskite privačius duomenis tik turėdami duomenų savininko leidimą. Vietiniams
modeliams šie apribojimai netaikomi — duomenys lieka kompiuteryje.

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
