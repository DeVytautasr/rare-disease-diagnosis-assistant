# Diegimas

Įrankis skaito struktūrinių variantų kandidatų rinkinius ir parodo, kokie
duomenys juos pagrindžia. Šis dokumentas skirtas žmogui, kuris šio kodo nerašė.

---

## Ko NEREIKIA

Prieš pradedant verta pasakyti, ko **nereikia** — tai dažniausiai ir atbaido:

| Nereikia | Kodėl |
|---|---|
| **Referencinio genomo** (FASTA, .fai, .dict) | Įrankis niekada neatidaro FASTA failo. BAM antraštė pati aprašo chromosomas. |
| **samtools ar bcftools** | Tą darbą atlieka Python biblioteka `pysam` tame pačiame procese. |
| **delly** | delly *pagamina* kandidatų rinkinį. Šis įrankis jį tik *skaito*. Žr. NAUDOJIMAS.md. |
| **Vaizdo plokštės (GPU)** | Nė vienai pagrindinei funkcijai. |
| **Interneto ryšio veikimo metu** | **MINIMAL lygiui — ne.** Genų pavadinimų paieškai — taip (neprivaloma). **IGV paveikslėliams — taip, būtinai** (žr. žemiau). |
| **Administratoriaus teisių** | Išskyrus vieną atvejį — žr. „Klaida: nepavyko sukurti aplinkos“. |

Reikia tik: **Python 3.10 arba naujesnio** ir maždaug **300 MB** vietos diske
Python paketams.

---

## Pirmiausia: ar jūsų kompiuteris tinka

| Sistema | Ar veiks | Pastaba |
|---|---|---|
| **Linux** (taip pat WSL) | taip | patikrinta |
| **macOS** (Intel ir Apple Silicon) | taip | yra paruošti paketai |
| **Windows be WSL** | **ne** | |

**Windows be WSL neveiks.** Biblioteka `pysam`, be kurios įrankis neskaito BAM
failų, **Windows sistemai paruoštų paketų neturi** — jos autoriai jų neleidžia.
Bandymas diegti baigsis nesėkme bandant kompiliuoti iš pirminio kodo.

Jeigu turite Windows kompiuterį, yra du keliai:

1. **Įsidiegti WSL** (`wsl --install` PowerShell'e administratoriaus teisėmis),
   paskui viską daryti jame kaip Linux sistemoje. Tai atskiras darbas.
2. **Nediegti nieko** ir pažiūrėti, kaip įrankis veikia, per ekrano bendrinimą.
   Dažniausiai to visiškai pakanka.

Jei tikslas — tiesiog pamatyti, ką įrankis rodo, **antrasis kelias yra
greitesnis ir patikimesnis**. Diegimas prasmingas tada, kai norite dirbti su
savo duomenimis.

---

## Diegimas

```bash
cd <katalogas-su-kodu>
bash install.sh
```

Skriptas:

1. patikrina Python versiją,
2. sukuria atskirą Python aplinką kataloge `.venv`,
3. įdiegia tris paketus (`pysam`, `fastmcp`, `requests`),
4. **pats pasitikrina** ir parodo, kas veiks, o kas ne.

## Kaip patikrinti, ar pavyko

Viena komanda:

```bash
.venv/bin/python -m stage1_igv_assistant.ui --check
```

Ji nieko nepaleidžia — tik parodo būseną ir baigiasi. Jei matote
`CHECK: MINIMAL tier OK`, įdiegta teisingai.

Pranešimas atrodo maždaug taip:

```
  available tier: MINIMAL
    [x] MINIMAL   filter chain, four evidence layers, hand-entered coordinates
    [ ] FULL      IGV panels — NOT FOUND; panels will report the failure instead of rendering
    [ ] COMPLETE  ollama not reachable; the chat panel will be unavailable
    [x] exclude template  /home/.../human.hg38.excl.tsv  [built-in default]
    [x] data directory    /home/.../public_data  [config file ...]
```

Laužtiniuose skliaustuose `[x]` reiškia „veikia“, `[ ]` — „nėra“.
**Tuščias langelis nėra klaida.** Trys lygiai yra savarankiški:

| Lygis | Ką gaunate | Ko reikia papildomai |
|---|---|---|
| **MINIMAL** | filtrų grandinė, keturi įrodymų sluoksniai, ranka įvesta koordinatė | nieko |
| **FULL** | tas pats + IGV paveikslėliai | IGV ir Java |
| **COMPLETE** | tas pats + pokalbis su vietiniu modeliu | `ollama` |

Jei turite tik MINIMAL — įrankis veikia visas, tik be paveikslėlių ir be pokalbio.

## Paleidimas

```bash
.venv/bin/python -m stage1_igv_assistant.ui
```

Naršyklėje atsidarykite **http://127.0.0.1:8765**. Serveris klausosi tik
šio kompiuterio — iš tinklo prie jo prieiti negalima.

Sustabdyti: `Ctrl+C`.

---

## Konfigūracija

Įrankiui reikia pasakyti, **kur yra duomenys**. Trys būdai, nuo paprasčiausio:

### 1. Nurodyti failus komandinėje eilutėje

```bash
.venv/bin/python -m stage1_igv_assistant.ui \
    --dataset MEGINYS=/kelias/iki/meginys.bam \
    --candidates MEGINYS=/kelias/iki/meginys.bcf
```

`MEGINYS` — jūsų pasirinkta **etiketė**. Naršyklė ir pokalbio modelis mato
**tik etiketę**, niekada failo kelio.

### 2. Konfigūracijos failas

```bash
cp sv-assistant.conf.example sv-assistant.conf
```

Atverkite ir pataisykite:

```ini
[paths]
data_dir = /kelias/iki/duomenu

[datasets]
MEGINYS = /kelias/iki/meginys.bam

[candidates]
MEGINYS = /kelias/iki/meginys.bcf
```

Santykiniai keliai skaičiuojami **nuo paties konfigūracijos failo katalogo**,
todėl demonstracinį rinkinį galima išskleisti bet kur.

Failas ieškomas eilės tvarka: `$SV_CONFIG`, `./sv-assistant.conf`,
`~/.config/sv-assistant/config.ini`.

### 3. Aplinkos kintamieji

| Kintamasis | Ką nurodo |
|---|---|
| `SV_DATA_DIR` | katalogas, kuriame ieškoti `.bam` ir `.bcf` |
| `SV_EXCLUDE_TEMPLATE` | delly „neįtraukiamų sričių“ šablonas |
| `IGV_PATH` | kelias iki `igv.sh` |
| `SV_OLLAMA_URL` | kur veikia `ollama` (numatyta `http://127.0.0.1:11434`) |
| `ANTHROPIC_API_KEY` | raktas Claude pokalbiui (neprivaloma) |

Pirmenybė: komandinė eilutė → aplinkos kintamasis → konfigūracijos failas →
numatytoji reikšmė.

---

## Demonstracinis rinkinys

Kataloge `demo_bundle` (**1,6 MB**) yra viskas, ko reikia pamatyti, kaip veikia:

```bash
SV_CONFIG=/kelias/iki/demo_bundle/sv-assistant.conf \
  .venv/bin/python -m stage1_igv_assistant.ui
```

Kas jame yra ir ko nėra — žr. NAUDOJIMAS.md skyrių „Demonstracinio rinkinio ribos“.

---

## Ar vadovei verta diegti?

Trumpas atsakymas: **jei reikia tik pamatyti, kaip veikia — nediekite.**
Parodykite per ekrano dalijimąsi. Užtrunka apie 10 minučių, jai nereikia
nieko įsirašyti ir nieko negali nepavykti.

Diegti verta **tik tada, jei ji norės pati dirbti su duomenimis** — pati
spaudyti filtrus, įvesti koordinates, grįžti prie to po savaitės.

### Ko jai NEREIKIA

Nereikia **referencinio genomo**, **delly**, **samtools**, **vaizdo plokštės**
ir **interneto ryšio darbo metu**. Nereikia ir patiems duomenims siųsti — visa
demonstracija telpa į 1,6 MB.

### Jei vis dėlto diegiama — trumpiausias kelias

Reikia tik dviejų dalykų: **Python 3.10 ar naujesnio** ir **interneto diegimo
metu** (tik tam, kad parsisiųstų Python paketus).

```bash
# 1. įdiegti
bash install.sh

# 2. patikrinti, ar pavyko — VIENA komanda, kuri viską pasako
.venv/bin/python -m stage1_igv_assistant.ui --check
```

Jei antroji komanda baigiasi eilute `CHECK: MINIMAL tier OK`, viskas gerai.
Jei ne — ta pati eilutė pasako, ko trūksta.

```bash
# 3. paleisti su demonstraciniu rinkiniu
SV_CONFIG=demo_bundle/sv-assistant.conf .venv/bin/python -m stage1_igv_assistant.ui
```

Tada naršyklėje atverti <http://127.0.0.1:8765>.

> **Windows:** įrankis rašytas Linux aplinkai. Windows kompiuteryje pirma reikia
> įjungti WSL (`wsl --install` administratoriaus teisėmis, po to perkrauti) ir
> visas komandas rašyti WSL lange. Jei tai skamba kaip kliūtis — grįžkite prie
> pirmos pastraipos: ekrano dalijimasis yra geresnis pasirinkimas.

---

## Kai kas nors nepavyksta

### `python3: command not found`
Python neįdiegtas. Ubuntu/WSL: `sudo apt install python3`.

### `ERROR: Python 3.10+ required`
Per sena versija. Įdiekite naujesnę ir nurodykite ją:
`PYTHON=/usr/bin/python3.12 bash install.sh`

### `ERROR: could not create a virtual environment`
Ubuntu ir WSL `venv` platinamas atskirai. Skriptas parašo tikslų paketo
pavadinimą, pvz.:

```bash
sudo apt install python3.12-venv
```

Tai **vienintelis** žingsnis, kuriam reikia administratoriaus teisių.
Paskui paleiskite `bash install.sh` iš naujo.

### `.venv` yra, bet diegimas nutrūksta
Nepavykęs ankstesnis bandymas paliko nebaigtą aplinką:

```bash
rm -rf .venv && bash install.sh
```

### `pysam` diegimas nepavyksta (klaida kompiliuojant)
Jūsų Python versijai dar nėra paruošto paketo. Naudokite Python 3.11 arba 3.12.

### Sąrašai naršyklėje tušti
Įrankis nerado duomenų. Patikrinkite pranešimą prie `data directory` — ten
matyti, kur buvo ieškota ir iš kur tas kelias paimtas. Žr. „Konfigūracija“.

### `[ ] exclude template`
Nėra delly neįtraukiamų sričių šablono. Įrankis veikia, tik filtrų grandinėje
tas žingsnis **nebus atliktas**, ir grandinė apie tai aiškiai parašys.
Tai nėra klaida — tiesiog vienu filtru mažiau.

### IGV paveikslėliai nesusikuria, nors IGV rastas

IGV pats **parsisiunčia genomo seką iš interneto** (`igv.org`) kiekvieną kartą,
jei ji nėra išsaugota vietoje. Jei to serverio pasiekti nepavyksta, IGV pakimba
ir įrankis po 120 s praneša `IGV timed out after 120s`.

Patikrinti:

```bash
curl -sS -o /dev/null -w '%{http_code}\n' https://igv.org/genomes/genomes.json
```

Jei atsakymo nėra — tai **ne šio įrankio klaida**. Visi skaičiai ir visi keturi
įrodymų sluoksniai veikia toliau; trūksta tik paveikslėlių.

Ilgalaikis sprendimas — vieną kartą paleisti IGV rankomis, pasirinkti `hg38`
genomą ir leisti jam išsisaugoti kataloge `~/igv/genomes`. Po to paveikslėliai
veiks ir be interneto.

IGV taip pat perspėja `IGV requires Java 17`, jei sistemoje yra naujesnė Java.
Iki šiol tai veikė ir su naujesne, bet jei paveikslėliai nesikuria — verta
pabandyti Java 17.

### `[ ] FULL — IGV not found`
Nėra IGV. Visa kita veikia; vietoj paveikslėlio matysite paaiškinimą, kur buvo
ieškota. **Skaičiai nuo to nepriklauso** — juos duoda skaičiavimo įrankiai, ne IGV.

### `[ ] COMPLETE — ollama not reachable`
Neveikia `ollama`. Pokalbio skydelio nebus; visa kita veikia.

### Prievadas 8765 užimtas
`.venv/bin/python -m stage1_igv_assistant.ui --port 8790`
