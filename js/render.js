import * as param from "./const.js"

const fixDiv = document.querySelector('.app__fix');
const tbody = document.querySelector('#data');

fixDiv.innerHTML = `
    <div><b>t1</b>${param.t1} [c]</div>
    <div><b>t2</b>${param.t2} [c]</div>
    <div><b>thet_torch</b>${param.thet_torch} [рад]</div>
    <div><b>thet_2</b>${param.thet_2} [рад]</div>
`
const tr_iter = (t,el, sign = 4) => {
    return `
        <td class="app__ceil">${t.toFixed(6)}</td>
        <td class="app__ceil">${el[0].toFixed(sign)}</td>
        <td class="app__ceil">${el[1].toFixed(sign)}</td>
        <td class="app__ceil">${el[2].toFixed(sign)}</td>
        <td class="app__ceil">${el[3].toFixed(sign)}</td>
        <td class="app__ceil">${el[4].toFixed(sign)}</td>
        <td class="app__ceil">${el[5].toFixed(sign)}</td>
        <td class="app__ceil">${el[6].toFixed(sign)}</td>
        <td class="app__ceil">${el[7].toFixed(sign)}</td>
        <td class="app__ceil">${el[8].toFixed(sign)}</td>
        <td class="app__ceil">${el[9].toFixed(sign)}</td>
        <td class="app__ceil">${el[10].toFixed(sign)}</td>
        <td class="app__ceil">${el[11].toFixed(sign)}</td>
    `
}

const html = (t, el) => {
    // const tr = document.createElement('tr')
    // tr.insertAdjacentHTML('beforeend', tr_iter(t, el))
    // tbody.insertAdjacentElement('beforeend', tr)
    console.log(`${t.toFixed(7)} ${el[0].toFixed(4)} ${el[1].toFixed(4)} ${el[2].toFixed(4)} ${el[3].toFixed(4)} ${el[4].toFixed(4)} ${el[5].toFixed(4)} ${el[6].toFixed(4)} ${el[7].toFixed(4)} ${el[8].toFixed(4)} ${el[9].toFixed(4)} ${el[10].toFixed(4)} ${el[11].toFixed(4)}`);
}

export {
    tr_iter,
    html
}